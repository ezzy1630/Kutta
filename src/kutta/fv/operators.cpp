#include "kutta/fv/operators.hpp"

#include <algorithm>
#include <array>
#include <cmath>

#include "kernels.hpp"
#include "kutta/mesh/detail.hpp"
#include "kutta/mesh/topology.hpp"

namespace kutta::fv {
using namespace detail;
namespace {
// One component at a time: three float copies would add half again to the
// Vec3d array's host footprint during setup.
DeviceVec3 vectors(Device& dev, const std::vector<Vec3d>& values) {
  std::vector<float> c(values.size());
  DeviceVec3 out;
  for (std::size_t i = 0; i < c.size(); ++i) c[i] = float(values[i].x);
  out.x = dev.upload<float>(c);
  for (std::size_t i = 0; i < c.size(); ++i) c[i] = float(values[i].y);
  out.y = dev.upload<float>(c);
  for (std::size_t i = 0; i < c.size(); ++i) c[i] = float(values[i].z);
  out.z = dev.upload<float>(c);
  return out;
}
// The rows of boundary faces, in boundary-face order, from per-entry rows.
// A boundary face has exactly one cell-face entry (its owner's).
std::vector<Vec3d> boundary_rows(const PolyMesh& m,const kutta::CellFaces& cf,const std::vector<Vec3d>& rows) {
  std::vector<Vec3d> result(m.n_faces()-m.n_internal_faces());
  for(Index j=0;j<Index(cf.faces.size());++j)
    if(cf.faces[j]>=m.n_internal_faces()) result[cf.faces[j]-m.n_internal_faces()]=rows[j];
  return result;
}
}

DeviceVec3 vector_field(Device& dev, std::size_t n, float value) {
  DeviceVec3 v{dev.alloc<float>(n),dev.alloc<float>(n),dev.alloc<float>(n)};
  dev.fill(v.x,value); dev.fill(v.y,value); dev.fill(v.z,value);
  return v;
}

Stencil make_stencil(Device& dev, const PolyMesh& m, const MeshGeometry& g) {
  Stencil s;
  const auto cf=build_cell_faces(m);
  const auto cn=build_cell_neighbours(m);
  auto& p=s.pattern;
  p.n_rows=m.n_cells;
  std::vector<Index> adjacent,entry;
  std::vector<Vec3d> delta,rf;
  std::vector<float> inverse(6*m.n_cells), diffusion;
  std::vector<Vec3d> nonorthogonal;
  for (Index c=0;c<m.n_cells;++c) {
    std::vector<Index> cols(cn.cells.begin()+cn.offsets[c],cn.cells.begin()+cn.offsets[c+1]);
    cols.insert(std::lower_bound(cols.begin(),cols.end(),c),c);
    for (auto j:cols) { if(j==c) p.diag.push_back(Index(p.col.size())); p.col.push_back(j); }
    p.row_ptr.push_back(Index(p.col.size()));
    std::array<double,6> a{};
    for (auto j=cf.offsets[c];j<cf.offsets[c+1];++j) {
      const Index f=cf.faces[j];
      const Index n=f<m.n_internal_faces() ? (m.owner[f]==c ? m.neighbour[f] : m.owner[f]) : -(f-m.n_internal_faces()+1);
      adjacent.push_back(n);
      entry.push_back(n<0 ? -1 : Index(std::lower_bound(p.col.begin()+p.row_ptr[c],p.col.end(),n)-p.col.begin()));
      const auto d=(n<0 ? g.face_centroid[f] : g.cell_centroid[n])-g.cell_centroid[c];
      const double w=1/dot(d,d);
      delta.push_back(d*w); rf.push_back(g.face_centroid[f]-g.cell_centroid[c]);
      a[0]+=w*d.x*d.x; a[1]+=w*d.x*d.y; a[2]+=w*d.x*d.z;
      a[3]+=w*d.y*d.y; a[4]+=w*d.y*d.z; a[5]+=w*d.z*d.z;
    }
    const double det=a[0]*(a[3]*a[5]-a[4]*a[4])-a[1]*(a[1]*a[5]-a[4]*a[2])+a[2]*(a[1]*a[4]-a[3]*a[2]);
    KUTTA_REQUIRE(det>1e-14,"FV LSQ stencil is rank deficient");
    const std::array<double,6> inv{(a[3]*a[5]-a[4]*a[4])/det,(a[2]*a[4]-a[1]*a[5])/det,
      (a[1]*a[4]-a[2]*a[3])/det,(a[0]*a[5]-a[2]*a[2])/det,
      (a[1]*a[2]-a[0]*a[4])/det,(a[0]*a[3]-a[1]*a[1])/det};
    for(int k=0;k<6;++k) inverse[c+k*m.n_cells]=float(inv[k]);
  }
  for(Index f=0;f<m.n_faces();++f) {
    const auto d=(f<m.n_internal_faces() ? g.cell_centroid[m.neighbour[f]] : g.face_centroid[f])-g.cell_centroid[m.owner[f]];
    const auto area=g.face_area[f];
    const double projection=dot(area,d);
    KUTTA_REQUIRE(projection>0,"FV requires positive face/cell projection");
    const double coeff=dot(area,area)/projection;
    diffusion.push_back(float(coeff)); nonorthogonal.push_back(area-d*coeff);
  }
  s.adjacent=dev.upload<Index>(adjacent); s.matrix_entry=dev.upload<Index>(entry);
  s.lsq_delta=vectors(dev,delta); s.face_delta=vectors(dev,rf);
  s.boundary_lsq_delta=vectors(dev,boundary_rows(m,cf,delta));
  s.lsq_inverse=dev.upload<float>(inverse); s.diffusion=dev.upload<float>(diffusion);
  s.nonorthogonal=vectors(dev,nonorthogonal);
  return s;
}

namespace {
// One cell's fit: its weighted rows d/|d|^2 appended to delta (zero for
// omitted boundary faces, in cell-face order) and the inverse of its normal
// matrix sum d d^T/|d|^2.
void cell_fit(const PolyMesh& m,const MeshGeometry& g,const kutta::CellFaces& cf,
    std::span<const BoundaryObservation> observations,Index c,std::vector<Vec3d>& delta,double inv[3][3]) {
  double a[3][3]{},q[3][3]{{1,0,0},{0,1,0},{0,0,1}};
  for(Index j=cf.offsets[c];j<cf.offsets[c+1];++j) {
    const auto f=cf.faces[j];
    const bool internal=f<m.n_internal_faces();
    if(!internal && observations[f-m.n_internal_faces()]==BoundaryObservation::Omitted) {delta.push_back({});continue;}
    const auto n=internal?(m.owner[f]==c?m.neighbour[f]:m.owner[f]):-1;
    auto d=(internal?g.cell_centroid[n]:g.face_centroid[f])-g.cell_centroid[c];
    if(!internal && observations[f-m.n_internal_faces()]==BoundaryObservation::NormalFootValue) {
      const auto area=g.face_area[f];d=area*(dot(d,area)/dot(area,area));
    }
    const double weight=1/dot(d,d),x[]{d.x,d.y,d.z};
    delta.push_back(d*weight);
    for(int r=0;r<3;++r) for(int k=0;k<3;++k) a[r][k]+=weight*x[r]*x[k];
  }
  // Symmetric Jacobi eigensolve. An inverse with unit eigenvalues only in
  // deficient directions is equivalent to adding orthogonal zero rows there;
  // unlike a diagonal ridge it does not bias any resolved gradient direction.
  for(int sweep=0;sweep<24;++sweep) {
    int p=0,r=1;
    for(int i=0;i<3;++i) for(int j=i+1;j<3;++j)
      if(std::abs(a[i][j])>std::abs(a[p][r])) {p=i;r=j;}
    if(std::abs(a[p][r])<1e-14) break;
    const double angle=.5*std::atan2(2*a[p][r],a[r][r]-a[p][p]);
    const double cs=std::cos(angle),sn=std::sin(angle);
    const double app=a[p][p],arr=a[r][r],apr=a[p][r];
    a[p][p]=cs*cs*app-2*cs*sn*apr+sn*sn*arr;
    a[r][r]=sn*sn*app+2*cs*sn*apr+cs*cs*arr;a[p][r]=a[r][p]=0;
    for(int i=0;i<3;++i) {
      if(i!=p && i!=r) {
        const double x=a[i][p],y=a[i][r];
        a[i][p]=a[p][i]=cs*x-sn*y;a[i][r]=a[r][i]=sn*x+cs*y;
      }
      const double x=q[i][p],y=q[i][r];q[i][p]=cs*x-sn*y;q[i][r]=sn*x+cs*y;
    }
  }
  for(int r=0;r<3;++r) for(int j=0;j<3;++j) inv[r][j]=0;
  for(int k=0;k<3;++k) {
    const double eigen=a[k][k]>1e-10?1/a[k][k]:1;
    for(int r=0;r<3;++r) for(int j=0;j<3;++j) inv[r][j]+=q[r][k]*eigen*q[j][k];
  }
}
}  // namespace

std::vector<double> boundary_extrapolation_gain(const PolyMesh& m,const MeshGeometry& g,
    std::span<const BoundaryObservation> observations) {
  KUTTA_REQUIRE(observations.size()==std::size_t(m.n_faces()-m.n_internal_faces()),
                "FV LSQ boundary mask size mismatch");
  const auto cf=build_cell_faces(m);
  const auto ni=m.n_internal_faces();
  std::vector<double> gain(m.n_faces()-ni);
  std::vector<Vec3d> delta;
  for(Index c=0;c<m.n_cells;++c) {
    double inv[3][3];
    delta.clear();cell_fit(m,g,cf,observations,c,delta,inv);
    for(Index j=cf.offsets[c];j<cf.offsets[c+1];++j) {
      const auto f=cf.faces[j];
      if(f<ni) continue;
      // r . M^-1 weighs each row's difference in the extrapolated face value.
      const auto r=g.face_centroid[f]-g.cell_centroid[c];
      const Vec3d w{r.x*inv[0][0]+r.y*inv[1][0]+r.z*inv[2][0],r.x*inv[0][1]+r.y*inv[1][1]+r.z*inv[2][1],
                    r.x*inv[0][2]+r.y*inv[1][2]+r.z*inv[2][2]};
      for(const auto& d:delta) gain[f-ni]+=std::abs(dot(w,d));
    }
  }
  return gain;
}

BoundaryLeastSquares boundary_least_squares(Device& dev,const PolyMesh& m,
    const MeshGeometry& g,std::span<const BoundaryObservation> observations,const Stencil* stencil) {
  KUTTA_REQUIRE(observations.size()==std::size_t(m.n_faces()-m.n_internal_faces()),
                "FV LSQ boundary mask size mismatch");
  const auto cf=build_cell_faces(m);
  std::vector<Vec3d> delta;
  std::vector<float> inverse(6*m.n_cells);
  for(Index c=0;c<m.n_cells;++c) {
    double inv[3][3];
    cell_fit(m,g,cf,observations,c,delta,inv);
    const double v[]{inv[0][0],inv[0][1],inv[0][2],inv[1][1],inv[1][2],inv[2][2]};
    for(int k=0;k<6;++k) inverse[c+k*m.n_cells]=float(v[k]);
  }
  BoundaryLeastSquares fit;
  fit.boundary_delta=vectors(dev,boundary_rows(m,cf,delta));
  fit.inverse=dev.upload<float>(inverse);
  if(stencil) {
    // Internal rows here and in make_stencil are the same double expression
    // (d/|d|^2 between centroids), so the stencil's floats are these rows.
    KUTTA_REQUIRE(stencil->lsq_delta.x.size()==delta.size(),"FV LSQ stencil does not match the mesh");
    fit.stencil_rows=&stencil->lsq_delta;
  } else fit.delta=vectors(dev,delta);
  return fit;
}

BoundaryLeastSquares green_gauss_cell(Device& dev,const PolyMesh& m,const MeshGeometry& g) {
  const auto cf=build_cell_faces(m);
  std::vector<Vec3d> delta;
  for(Index c=0;c<m.n_cells;++c) {
    const double volume=g.cell_volume[c];
    KUTTA_REQUIRE(volume>0,"FV Green-Gauss needs positive cell volumes");
    for(Index j=cf.offsets[c];j<cf.offsets[c+1];++j) {
      const auto f=cf.faces[j];
      const auto outward=g.face_area[f]*double(cf.signs[j]);
      delta.push_back(outward*((f<m.n_internal_faces()?.5:1.)/volume));
    }
  }
  std::vector<float> inverse(6*m.n_cells,0.f);
  for(Index c=0;c<m.n_cells;++c) inverse[c]=inverse[c+3*m.n_cells]=inverse[c+5*m.n_cells]=1.f;
  BoundaryLeastSquares fit;
  fit.delta=vectors(dev,delta);
  fit.boundary_delta=vectors(dev,boundary_rows(m,cf,delta));
  fit.inverse=dev.upload<float>(inverse);
  return fit;
}

FaceJumpFit face_jump_least_squares(Device& dev,const PolyMesh& m,const MeshGeometry& g) {
  const auto cf=build_cell_faces(m);std::vector<Vec3d> delta;
  std::vector<float> inverse(6*m.n_cells);
  for(Index c=0;c<m.n_cells;++c) {
    std::array<double,6> a{};
    for(Index j=cf.offsets[c];j<cf.offsets[c+1];++j) {
      const auto r=g.face_centroid[cf.faces[j]]-g.cell_centroid[c];const double w=1/dot(r,r);
      delta.push_back(r*w);
      a[0]+=w*r.x*r.x;a[1]+=w*r.x*r.y;a[2]+=w*r.x*r.z;
      a[3]+=w*r.y*r.y;a[4]+=w*r.y*r.z;a[5]+=w*r.z*r.z;
    }
    const double det=a[0]*(a[3]*a[5]-a[4]*a[4])-a[1]*(a[1]*a[5]-a[4]*a[2])+a[2]*(a[1]*a[4]-a[3]*a[2]);
    KUTTA_REQUIRE(det>1e-14,"FV face-jump fit is rank deficient");
    const std::array<double,6> inv{(a[3]*a[5]-a[4]*a[4])/det,(a[2]*a[4]-a[1]*a[5])/det,
      (a[1]*a[4]-a[2]*a[3])/det,(a[0]*a[5]-a[2]*a[2])/det,
      (a[1]*a[2]-a[0]*a[4])/det,(a[0]*a[3]-a[1]*a[1])/det};
    for(int k=0;k<6;++k) inverse[c+k*m.n_cells]=float(inv[k]);
  }
  return {vectors(dev,delta),dev.upload<float>(inverse)};
}
void jump_corrected_gradient(Device& dev,const DeviceMesh& m,const Stencil& s,
    const FaceJumpFit& fit,const Buffer<float>& cell,const Buffer<float>& bc,
    const DeviceVec3& base,DeviceVec3& out) {
  kern::fv_jump_gradient(dev,m.n_cells,m.cell_face_offsets,m.cell_faces,s.adjacent,m.cell_face_sign,
    fit.delta.x,fit.delta.y,fit.delta.z,fit.inverse,s.face_delta.x,s.face_delta.y,s.face_delta.z,
    m.owner_to_face.x,m.owner_to_face.y,m.owner_to_face.z,
    m.neighbour_to_face.x,m.neighbour_to_face.y,m.neighbour_to_face.z,
    cell,bc,base.x,base.y,base.z,out.x,out.y,out.z,m.n_cells);
}

void gradient(Device& dev,const DeviceMesh& m,const Stencil& s,
    const BoundaryLeastSquares& lsq,const Buffer<float>& cell,
    const Buffer<float>& bc,DeviceVec3& out) {
  const auto& rows=lsq.rows();
  kern::fv_gradient(dev,m.n_cells,m.cell_face_offsets,s.adjacent,rows.x,rows.y,rows.z,
    lsq.boundary_delta.x,lsq.boundary_delta.y,lsq.boundary_delta.z,lsq.inverse,cell,bc,out.x,out.y,out.z,m.n_cells);
}
void vector_gradient(Device& dev,const DeviceMesh& m,const Stencil& s,
    const BoundaryLeastSquares& lsq,const DeviceVec3& cell,const DeviceVec3& bc,
    DeviceVec3& u,DeviceVec3& v,DeviceVec3& w) {
  const auto& rows=lsq.rows();
  kern::fv_vector_gradient(dev,m.n_cells,m.cell_face_offsets,s.adjacent,rows.x,rows.y,rows.z,
    lsq.boundary_delta.x,lsq.boundary_delta.y,lsq.boundary_delta.z,lsq.inverse,
    cell.x,cell.y,cell.z,bc.x,bc.y,bc.z,u.x,u.y,u.z,v.x,v.y,v.z,w.x,w.y,w.z,m.n_cells);
}

void gradient(Device& dev, const DeviceMesh& m, const Stencil& s,
              const Buffer<float>& cell, const Buffer<float>& boundary, DeviceVec3& out) {
  KUTTA_REQUIRE(cell.size()==std::size_t(m.n_cells) && boundary.size()==std::size_t(m.n_boundary_faces()) &&
    out.x.size()==cell.size() && out.y.size()==cell.size() && out.z.size()==cell.size(),"FV gradient field size mismatch");
  kern::fv_gradient(dev,m.n_cells,m.cell_face_offsets,s.adjacent,s.lsq_delta.x,s.lsq_delta.y,s.lsq_delta.z,
    s.boundary_lsq_delta.x,s.boundary_lsq_delta.y,s.boundary_lsq_delta.z,
    s.lsq_inverse,cell,boundary,out.x,out.y,out.z,m.n_cells);
}
// With M = max|phi|: max(phi) = max|phi + M| - M and min(phi) = M - max|M - phi|.
std::vector<FieldRange> field_ranges(Device& dev, std::span<const Buffer<float>* const> fields, std::size_t n,
                                     std::span<Buffer<float>> scratch) {
  KUTTA_REQUIRE(scratch.size()>=2*fields.size(),"FV field range scratch mismatch");
  std::vector<FieldRange> out(fields.size());
  if(n==0 || fields.empty()) return out;
  std::vector<Device::HostReduction> r;
  for(std::size_t k=0;k<fields.size();++k) {
    KUTTA_REQUIRE(n<=fields[k]->size() && n<=scratch[2*k].size() && n<=scratch[2*k+1].size(),"FV field range size mismatch");
    r.push_back(Device::HostReduction::max_abs(*fields[k],n));
  }
  const auto bound=dev.reduce_many(r);
  r.clear();
  for(std::size_t k=0;k<fields.size();++k) {
    kern::fv_range_shift(dev,n,*fields[k],scratch[2*k],scratch[2*k+1],float(bound[k]));
    r.push_back(Device::HostReduction::max_abs(scratch[2*k],n));r.push_back(Device::HostReduction::max_abs(scratch[2*k+1],n));
  }
  const auto m=dev.reduce_many(r);
  for(std::size_t k=0;k<fields.size();++k) {
    const auto b=float(bound[k]);out[k]={std::max(0.0,m[2*k]+m[2*k+1]-2.0*double(b)),b};
  }
  return out;
}

void limit_gradient(Device& dev, const DeviceMesh& m, const Stencil& s,
                    const Buffer<float>& cell, const Buffer<float>& boundary,
                    const DeviceVec3& g, Buffer<float>& limiter,
                    const Buffer<float>* low, const Buffer<Index>* boundary_kind, float epsilon2) {
  KUTTA_REQUIRE(limiter.size()==cell.size() && (!low || (low->size()==cell.size() && boundary_kind &&
    boundary_kind->size()==boundary.size())),"FV limiter size mismatch");
  if(low) kern::fv_paired_limiter(dev,m.n_cells,m.cell_face_offsets,s.adjacent,s.face_delta.x,s.face_delta.y,
    s.face_delta.z,cell,*low,boundary,*boundary_kind,g.x,g.y,g.z,limiter,epsilon2);
  else kern::fv_limiter(dev,m.n_cells,m.cell_face_offsets,s.adjacent,s.face_delta.x,s.face_delta.y,s.face_delta.z,
    cell,boundary,g.x,g.y,g.z,limiter,epsilon2);
}

void vector_gradient(Device& dev,const DeviceMesh& m,const Stencil& s,
    const DeviceVec3& cell,const DeviceVec3& bc,DeviceVec3& u,DeviceVec3& v,DeviceVec3& w) {
  bool sizes=true;
  const Buffer<float>* fields[]{&cell.x,&cell.y,&cell.z,&u.x,&u.y,&u.z,&v.x,&v.y,&v.z,&w.x,&w.y,&w.z};
  for(const auto* field:fields)
    sizes=sizes && field->size()==std::size_t(m.n_cells);
  for(const auto* field:{&bc.x,&bc.y,&bc.z}) sizes=sizes && field->size()==std::size_t(m.n_boundary_faces());
  KUTTA_REQUIRE(sizes,"FV vector gradient size mismatch");
  kern::fv_vector_gradient(dev,m.n_cells,m.cell_face_offsets,s.adjacent,
    s.lsq_delta.x,s.lsq_delta.y,s.lsq_delta.z,
    s.boundary_lsq_delta.x,s.boundary_lsq_delta.y,s.boundary_lsq_delta.z,s.lsq_inverse,
    cell.x,cell.y,cell.z,bc.x,bc.y,bc.z,u.x,u.y,u.z,v.x,v.y,v.z,w.x,w.y,w.z,m.n_cells);
}
}  // namespace kutta::fv
