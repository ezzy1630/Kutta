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
DeviceVec3 vectors(Device& dev, const std::vector<Vec3d>& values) {
  std::vector<float> x, y, z;
  for (auto v : values) { x.push_back(float(v.x)); y.push_back(float(v.y)); z.push_back(float(v.z)); }
  return {dev.upload<float>(x),dev.upload<float>(y),dev.upload<float>(z)};
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
  s.lsq_inverse=dev.upload<float>(inverse); s.diffusion=dev.upload<float>(diffusion);
  s.nonorthogonal=vectors(dev,nonorthogonal);
  return s;
}

void gradient(Device& dev, const DeviceMesh& m, const Stencil& s,
              const Buffer<float>& cell, const Buffer<float>& boundary, DeviceVec3& out) {
  KUTTA_REQUIRE(cell.size()==std::size_t(m.n_cells) && boundary.size()==std::size_t(m.n_boundary_faces()) &&
    out.x.size()==cell.size() && out.y.size()==cell.size() && out.z.size()==cell.size(),"FV gradient field size mismatch");
  kern::fv_gradient(dev,m.n_cells,m.cell_face_offsets,s.adjacent,s.lsq_delta.x,s.lsq_delta.y,s.lsq_delta.z,
    s.lsq_inverse,cell,boundary,out.x,out.y,out.z,m.n_cells);
}
void limit_gradient(Device& dev, const DeviceMesh& m, const Stencil& s,
                    const Buffer<float>& cell, const Buffer<float>& boundary,
                    const DeviceVec3& g, Buffer<float>& limiter) {
  KUTTA_REQUIRE(limiter.size()==cell.size(),"FV limiter size mismatch");
  kern::fv_limiter(dev,m.n_cells,m.cell_face_offsets,s.adjacent,s.face_delta.x,s.face_delta.y,s.face_delta.z,
    cell,boundary,g.x,g.y,g.z,limiter);
}
}  // namespace kutta::fv
