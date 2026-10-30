#include "kutta/fv/face_projection.hpp"
#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <map>
#include <set>
#include <limits>
#include <thread>
#include "kernels.hpp"
#include "kutta/mesh/detail.hpp"
#include "kutta/mesh/topology.hpp"

namespace kutta::fv {
using namespace detail;
namespace {
constexpr int terms=10;
using Row=std::array<double,terms>;
using Matrix=std::array<Row,terms>;
Row basis(Vec3d r) {return {r.x,r.y,r.z,.5*r.x*r.x,r.x*r.y,r.x*r.z,.5*r.y*r.y,r.y*r.z,.5*r.z*r.z,1};}
Matrix inverse(Matrix a) {
  Matrix q{};for(int i=0;i<terms;++i) q[i][i]=1;
  for(int iteration=0;iteration<terms*terms*16;++iteration) {
    int p=0,r=1;double scale=0;
    for(int i=0;i<terms;++i) {
      scale=std::max(scale,std::abs(a[i][i]));
      for(int j=i+1;j<terms;++j) if(std::abs(a[i][j])>std::abs(a[p][r])) {p=i;r=j;}
    }
    if(std::abs(a[p][r])<1e-13*scale) break;
    const double angle=.5*std::atan2(2*a[p][r],a[r][r]-a[p][p]);
    const double c=std::cos(angle),s=std::sin(angle),app=a[p][p],arr=a[r][r],apr=a[p][r];
    a[p][p]=c*c*app-2*c*s*apr+s*s*arr;a[r][r]=s*s*app+2*c*s*apr+c*c*arr;a[p][r]=a[r][p]=0;
    for(int i=0;i<terms;++i) {
      if(i!=p && i!=r) {const double x=a[i][p],y=a[i][r];a[i][p]=a[p][i]=c*x-s*y;a[i][r]=a[r][i]=s*x+c*y;}
      const double x=q[i][p],y=q[i][r];q[i][p]=c*x-s*y;q[i][r]=s*x+c*y;
    }
  }
  double largest=0;for(int i=0;i<terms;++i) largest=std::max(largest,a[i][i]);
  Matrix out{};
  for(int k=0;k<terms;++k) if(a[k][k]>largest*1e-11)
    for(int i=0;i<terms;++i) for(int j=0;j<terms;++j) out[i][j]+=q[i][k]*q[j][k]/a[k][k];
  return out;
}
Row product(const Matrix& a,const Row& b) {
  Row out{};for(int i=0;i<terms;++i) for(int j=0;j<terms;++j) out[i]+=a[i][j]*b[j];return out;
}
double dot_row(const Row& a,const Row& b) {double out=0;for(int i=0;i<terms;++i) out+=a[i]*b[i];return out;}
}
ProjectedFaceReconstruction::ProjectedFaceReconstruction(Device& dev,const PolyMesh& mesh,const MeshGeometry& g,
    std::span<const BoundaryObservation> observations,std::span<const Index> nonsmooth_cells):dev_(dev),faces_(mesh.n_faces()),internal_(mesh.n_internal_faces()) {
  const auto internal=mesh.n_internal_faces();
  KUTTA_REQUIRE(nonsmooth_cells.empty() || nonsmooth_cells.size()==std::size_t(mesh.n_cells),"face projection smoothness mask mismatch");
  KUTTA_REQUIRE(observations.size()==std::size_t(faces_-internal),"face projection boundary mask mismatch");
  const auto cn=build_cell_neighbours(mesh);const auto cf=build_cell_faces(mesh);
  std::vector<Index> offsets{0},points,mask(faces_,0);
  std::vector<float> wx,wy,wz,wflux,wvalue,wlink;
  for(Index face=0;face<faces_;++face) {
    if(face>=internal && observations[face-internal]!=BoundaryObservation::FaceValue) {offsets.push_back(Index(points.size()));continue;}
    const auto owner=mesh.owner[face];
    const auto face_area=g.face_area[face],unit_normal=face_area/norm(face_area);
    const auto rp=g.face_centroid[face]-g.cell_centroid[owner];
    const auto rn=face<internal?g.cell_centroid[mesh.neighbour[face]]-g.face_centroid[face]:Vec3d{};
    const double dp=dot(rp,unit_normal),dn=dot(rn,unit_normal);
    // The original stencil already is the normal projection at roundoff.
    // Include both cells in the coordinate scale so ownership cannot alter it.
    const double roundoff=256*std::numeric_limits<double>::epsilon()*std::max({1.,norm(g.face_centroid[face]),norm(g.cell_centroid[owner]),face<internal?norm(g.cell_centroid[mesh.neighbour[face]]):0.});
    if(norm(rp-unit_normal*dp)<=roundoff && norm(rn-unit_normal*dn)<=roundoff) {offsets.push_back(Index(points.size()));continue;}


    bool accepted=false;
    // Overlap gives the two adjacent cells the same smooth reconstruction.
    // Expand only when the requested functionals are rank deficient.
    for(int rings=3;rings<=5 && !accepted;++rings) {
      std::set<Index> cloud{owner};if(face<internal) cloud.insert(mesh.neighbour[face]);
      for(int depth=0;depth<rings;++depth) {
        const auto old=cloud;
        for(auto c:old) for(Index j=cn.offsets[c];j<cn.offsets[c+1];++j) cloud.insert(cn.cells[j]);
      }
      if(!nonsmooth_cells.empty() && std::any_of(cloud.begin(),cloud.end(),[&](Index c){return nonsmooth_cells[c]!=0;})) break;
      std::set<Index> samples=cloud;
      for(auto c:cloud) for(Index j=cf.offsets[c];j<cf.offsets[c+1];++j) {
        const auto f=cf.faces[j];if(f>=internal && observations[f-internal]==BoundaryObservation::FaceValue) samples.insert(-(f-internal+1));
      }
      if(samples.empty()) continue;
      auto point=[&](Index i) {return i<0?g.face_centroid[internal-i-1]:g.cell_centroid[i];};
      const auto fit_origin=g.face_centroid[face];
      Vec3d scale{};
      for(auto i:samples) {const auto r=point(i)-fit_origin;scale={std::max(scale.x,std::abs(r.x)),std::max(scale.y,std::abs(r.y)),std::max(scale.z,std::abs(r.z))};}
      const double length=std::max({scale.x,scale.y,scale.z});
      if(!(length>0) || !std::isfinite(length)) continue;
      const bool resolved[]{scale.x>=length*1e-10,scale.y>=length*1e-10,scale.z>=length*1e-10};
      if(scale.x<length*1e-10) scale.x=length;if(scale.y<length*1e-10) scale.y=length;if(scale.z<length*1e-10) scale.z=length;
      auto scaled=[&](Vec3d r) {return Vec3d{r.x/scale.x,r.y/scale.y,r.z/scale.z};};
      Matrix a{};std::vector<Row> rows;std::vector<double> weights;
      for(auto i:samples) {
        const auto r=scaled(point(i)-fit_origin);const auto b=basis(r);const auto sample=point(i),other=face<internal?g.cell_centroid[mesh.neighbour[face]]:g.face_centroid[face];
        const double w=1/(dot(sample-g.cell_centroid[owner],sample-g.cell_centroid[owner])+dot(sample-other,sample-other));
        rows.push_back(b);weights.push_back(w);
        for(int k=0;k<terms;++k) for(int j=0;j<terms;++j) a[k][j]+=w*b[k]*b[j];
      }
      const auto inv=inverse(a);const auto r=scaled(g.face_centroid[face]-fit_origin);
      const Row targets[]{Row{1,0,0,r.x,r.y,r.z,0,0,0},Row{0,1,0,0,r.x,0,r.y,r.z,0},Row{0,0,1,0,0,r.x,0,r.y,r.z}};
      const auto area=g.face_area[face];
      const auto link=(face<internal?g.cell_centroid[mesh.neighbour[face]]:g.face_centroid[face])-g.cell_centroid[owner];
      const double alpha=dot(area,area)/dot(area,link);
      const auto owner_row=basis(scaled(g.cell_centroid[owner]-fit_origin));
      auto jump=basis(scaled((face<internal?g.cell_centroid[mesh.neighbour[face]]:g.face_centroid[face])-fit_origin));
      for(int k=0;k<terms;++k) jump[k]-=owner_row[k];
      // Project P and N tangentially onto the normal through F. Preserve the
      // normal finite difference (including its curvature error); reconstruct
      // only the shifts q(P*)-q(P), q(N*)-q(N). It vanishes on orthogonal faces.
      auto pfoot=basis(scaled(g.face_centroid[face]-unit_normal*dp-fit_origin));
      auto nfoot=basis(scaled(g.face_centroid[face]+unit_normal*dn-fit_origin));
      for(int k=0;k<terms;++k) {pfoot[k]-=owner_row[k];nfoot[k]-=owner_row[k];}
      Row correction{};for(int k=0;k<terms;++k) correction[k]=alpha*(nfoot[k]-jump[k]-pfoot[k]);
      const auto recovered=product(a,product(inv,correction));double error=0,size=0;
      for(int k=0;k<terms;++k) {error=std::max(error,std::abs(recovered[k]-correction[k]));size=std::max(size,std::abs(correction[k]));}
      if(error>1e-7*size) continue;
      Row face_value{};
      for(int k=0;k<terms;++k) face_value[k]=(dn*pfoot[k]+dp*nfoot[k])/(dp+dn);
      const auto recovered_value=product(a,product(inv,face_value));
      double value_error=0,value_size=0;
      for(int k=0;k<terms;++k) {value_error=std::max(value_error,std::abs(recovered_value[k]-face_value[k]));value_size=std::max(value_size,std::abs(face_value[k]));}
      if(value_error>1e-7*value_size) continue;
      bool gradient_resolved=true;
      for(int axis=0;axis<3;++axis) if(resolved[axis]) {
        const auto recovered_gradient=product(a,product(inv,targets[axis]));double residual=0,magnitude=0;
        for(int k=0;k<terms;++k) {residual=std::max(residual,std::abs(recovered_gradient[k]-targets[axis][k]));magnitude=std::max(magnitude,std::abs(targets[axis][k]));}
        if(residual>1e-7*magnitude) gradient_resolved=false;
      }
      if(!gradient_resolved) continue;
      const auto x=product(inv,targets[0]),y=product(inv,targets[1]),z=product(inv,targets[2]);
      const auto flux=product(inv,correction);
      const auto value=product(inv,face_value),link_value=product(inv,jump);
      std::size_t j=0;
      for(auto i:samples) {
        points.push_back(i);wvalue.push_back(float(weights[j]*dot_row(value,rows[j])));
        {
          wx.push_back(float(weights[j]*dot_row(x,rows[j])/scale.x));wy.push_back(float(weights[j]*dot_row(y,rows[j])/scale.y));wz.push_back(float(weights[j]*dot_row(z,rows[j])/scale.z));
          wflux.push_back(float(weights[j]*dot_row(flux,rows[j])));wlink.push_back(float(weights[j]*dot_row(link_value,rows[j])));
        }
        ++j;
      }
      mask[face]=1;accepted=true;++corrected_;
    }
    if(!accepted) ++fallback_;
    offsets.push_back(Index(points.size()));
  }
  offsets_=dev.upload<Index>(offsets);points_=dev.upload<Index>(points);owner_=dev.upload<Index>(mesh.owner);mask_=dev.upload<Index>(mask);
  value_weights_=dev.upload<float>(wvalue);
  {
    weights_={dev.upload<float>(wx),dev.upload<float>(wy),dev.upload<float>(wz)};
    diffusion_weights_=dev.upload<float>(wflux);link_weights_=dev.upload<float>(wlink);
    neighbour_=dev.upload<Index>(mesh.neighbour);
  }
}
void ProjectedFaceReconstruction::diffusion_correction(const Buffer<float>& field,const Buffer<float>& boundary,Buffer<float>& out) const {
  kern::fv_projected_diffusion_correction(dev_,faces_,offsets_,points_,owner_,neighbour_,diffusion_weights_,field,boundary,out,internal_);
}
void ProjectedFaceReconstruction::momentum_values(const Buffer<float>& field,const Buffer<float>& boundary,
    const Buffer<float>& diagonal,Buffer<float>& out,const Buffer<float>* cell_limiter) const {
  kern::fv_projected_momentum_faces(dev_,internal_,offsets_,points_,owner_,neighbour_,mask_,value_weights_,link_weights_,field,boundary,diagonal,cell_limiter?*cell_limiter:diagonal,out,cell_limiter?1:0);
}
void ProjectedFaceReconstruction::gradient(const Buffer<float>& field,const Buffer<float>& boundary,DeviceVec3& out) const {
  kern::fv_projected_face_gradient(dev_,faces_,offsets_,points_,owner_,neighbour_,weights_.x,weights_.y,weights_.z,mask_,
    field,boundary,out.x,out.y,out.z,internal_);
}

void ProjectedFaceReconstruction::flux_gradient(const Buffer<float>& field,const Buffer<float>& boundary,
    const DeviceMesh& mesh,const Stencil& stencil,DeviceVec3& out) const {
  gradient(field,boundary,out);
  kern::fv_flux_consistent_gradient(dev_,faces_,offsets_,points_,owner_,neighbour_,mask_,diffusion_weights_,field,boundary,
    mesh.area.x,mesh.area.y,mesh.area.z,stencil.diffusion,out.x,out.y,out.z,internal_);
}
} // namespace kutta::fv

namespace kutta::fv {
using namespace detail;
namespace {
// Basis about the old origin, re-expressed about the new one; d=old-new.
Row shift(const Row& m,Vec3d d) {
  const double c=m[9];
  return {m[0]+d.x*c,m[1]+d.y*c,m[2]+d.z*c,m[3]+d.x*m[0]+.5*d.x*d.x*c,m[4]+d.y*m[0]+d.x*m[1]+d.x*d.y*c,
    m[5]+d.z*m[0]+d.x*m[2]+d.x*d.z*c,m[6]+d.y*m[1]+.5*d.y*d.y*c,m[7]+d.z*m[1]+d.y*m[2]+d.y*d.z*c,
    m[8]+d.z*m[2]+.5*d.z*d.z*c,c};
}
}
WalkFaceReconstruction::WalkFaceReconstruction(Device& dev,const PolyMesh& mesh,const MeshGeometry& g,
    const DeviceMesh& device_mesh,const Stencil& stencil,std::span<const BoundaryObservation> observations,
    std::span<const Index> nonsmooth_cells):dev_(dev),mesh_(device_mesh),stencil_(stencil) {
  const Index cells=mesh.n_cells,faces=mesh.n_faces(),internal=mesh.n_internal_faces();
  KUTTA_REQUIRE(device_mesh.n_cells==cells && device_mesh.n_faces==faces,"walk reconstruction mesh mismatch");
  KUTTA_REQUIRE(observations.size()==std::size_t(faces-internal),"walk reconstruction boundary mask mismatch");
  KUTTA_REQUIRE(nonsmooth_cells.empty() || nonsmooth_cells.size()==std::size_t(cells),"walk reconstruction smoothness mask mismatch");
  const auto cf=build_cell_faces(mesh);
  std::vector<Index> observed(faces-internal);
  for(Index b=0;b<faces-internal;++b) observed[b]=observations[b]==BoundaryObservation::FaceValue;
  const auto other=[&](Index c,Index f) {return mesh.owner[f]==c?mesh.neighbour[f]:mesh.owner[f];};
  // Geometry moments G^k (moments of the constant one), exactly as the device
  // hops accumulate them: one walk step per internal face of the cell.
  std::vector<std::vector<Row>> geometry(hops+1,std::vector<Row>(cells));
  for(Index c=0;c<cells;++c) {
    Row m{};m[9]=1;
    for(Index j=cf.offsets[c];j<cf.offsets[c+1];++j) {
      const auto f=cf.faces[j];if(f<internal || !observed[f-internal]) continue;
      const auto b=basis(g.face_centroid[f]-g.cell_centroid[c]);for(int k=0;k<terms;++k) m[k]+=b[k];
    }
    geometry[0][c]=m;
  }
  for(int level=0;level<hops;++level) for(Index c=0;c<cells;++c) {
    auto m=geometry[level][c];
    for(Index j=cf.offsets[c];j<cf.offsets[c+1];++j) {
      const auto f=cf.faces[j];if(f>=internal) continue;const auto d=other(c,f);
      const auto s=shift(geometry[level][d],g.cell_centroid[d]-g.cell_centroid[c]);for(int k=0;k<terms;++k) m[k]+=s[k];
    }
    geometry[level+1][c]=m;
  }
  // Each face's fit is independent; workers fill per-face results and the
  // slots are numbered serially afterwards, so threading never changes them.
  struct Fit {bool corrected=false,fallback=false;std::array<float,3*terms> c{};};
  std::vector<Fit> fits(faces);
  const auto fit_face=[&](Index face,std::map<Index,double>& walks,std::map<Index,double>& next) {
    Fit out;
    if(face>=internal && !observed[face-internal]) return out;
    const auto owner=mesh.owner[face];const bool interior=face<internal;
    const auto area=g.face_area[face],unit=area/norm(area),origin=g.face_centroid[face];
    const auto rp=origin-g.cell_centroid[owner],rn=interior?g.cell_centroid[mesh.neighbour[face]]-origin:Vec3d{};
    const double dp=dot(rp,unit),dn=dot(rn,unit);
    const double roundoff=256*std::numeric_limits<double>::epsilon()*std::max({1.,norm(origin),norm(g.cell_centroid[owner]),interior?norm(g.cell_centroid[mesh.neighbour[face]]):0.});
    if(norm(rp-unit*dp)<=roundoff && norm(rn-unit*dn)<=roundoff) return out;
    out.fallback=true;
    // Normal matrix sum_j w_j b_j b_j^T about the face centroid.
    Matrix a{};
    for(Index start:{owner,interior?mesh.neighbour[face]:Index(-1)}) {
      if(start<0) continue;
      walks.clear();walks[start]=1;
      for(int level=0;level<hops;++level) {
        next=walks;
        for(auto [c,x]:walks) for(Index j=cf.offsets[c];j<cf.offsets[c+1];++j) if(cf.faces[j]<internal) next[other(c,cf.faces[j])]+=x;
        walks.swap(next);
      }
      for(auto [c,x]:walks) {
        if(!nonsmooth_cells.empty() && nonsmooth_cells[c]!=0) return out;
        const auto add=[&](Vec3d point) {const auto b=basis(point-origin);for(int k=0;k<terms;++k) for(int l=0;l<terms;++l) a[k][l]+=x*b[k]*b[l];};
        add(g.cell_centroid[c]);
        for(Index j=cf.offsets[c];j<cf.offsets[c+1];++j) {const auto f=cf.faces[j];if(f>=internal && observed[f-internal]) add(g.face_centroid[f]);}
      }
    }
    if(!(a[9][9]>0)) return out;
    // Scale each axis by its weighted RMS sample offset for conditioning.
    std::array<double,3> scale{std::sqrt(a[0][0]/a[9][9]),std::sqrt(a[1][1]/a[9][9]),std::sqrt(a[2][2]/a[9][9])};
    const double length=*std::max_element(scale.begin(),scale.end());
    if(!(length>0) || !std::isfinite(length)) return out;
    for(auto& v:scale) if(v<length*1e-10) v=length;
    const Row e{scale[0],scale[1],scale[2],scale[0]*scale[0],scale[0]*scale[1],scale[0]*scale[2],scale[1]*scale[1],scale[1]*scale[2],scale[2]*scale[2],1};
    Matrix scaled{};for(int k=0;k<terms;++k) for(int l=0;l<terms;++l) scaled[k][l]=a[k][l]/(e[k]*e[l]);
    const auto inv=inverse(scaled);
    const auto other_point=interior?g.cell_centroid[mesh.neighbour[face]]:origin;
    const double alpha=dot(area,area)/dot(area,other_point-g.cell_centroid[owner]);
    const auto owner_row=basis(g.cell_centroid[owner]-origin);
    Row pfoot=basis(unit*(-dp)),nfoot=basis(unit*dn),jump=basis(other_point-origin);
    for(int k=0;k<terms;++k) {pfoot[k]-=owner_row[k];nfoot[k]-=owner_row[k];jump[k]-=owner_row[k];}
    Row targets[3]{};
    for(int k=0;k<terms;++k) {
      targets[0][k]=alpha*(nfoot[k]-jump[k]-pfoot[k]);
      targets[1][k]=(dn*pfoot[k]+dp*nfoot[k])/(dp+dn);
      targets[2][k]=jump[k];
    }
    for(int t=0;t<3;++t) {
      Row target{};for(int k=0;k<terms;++k) target[k]=targets[t][k]/e[k];
      const auto solution=product(inv,target),recovered=product(scaled,solution);
      double error=0,size=0;
      for(int k=0;k<terms;++k) {error=std::max(error,std::abs(recovered[k]-target[k]));size=std::max(size,std::abs(target[k]));}
      if(error>1e-7*size) return out;
      for(int k=0;k<terms;++k) out.c[t*terms+k]=float(solution[k]/e[k]);
    }
    out.corrected=true;out.fallback=false;return out;
  };
  {
    std::atomic<Index> cursor{0};constexpr Index chunk=256;
    const auto work=[&] {
      std::map<Index,double> walks,next;
      for(Index begin;(begin=cursor.fetch_add(chunk))<faces;)
        for(Index face=begin;face<std::min(faces,begin+chunk);++face) fits[face]=fit_face(face,walks,next);
    };
    std::vector<std::thread> workers(std::max(1u,std::min(std::thread::hardware_concurrency(),32u))-1);
    for(auto& w:workers) w=std::thread(work);
    work();for(auto& w:workers) w.join();
  }
  std::vector<Index> slot(faces,-1),mask(faces,0);
  std::vector<std::array<float,3*terms>> coefficients;
  for(Index face=0;face<faces;++face) {
    if(fits[face].fallback) ++fallback_;
    if(!fits[face].corrected) continue;
    slot[face]=Index(coefficients.size());mask[face]=1;coefficients.push_back(fits[face].c);++corrected_;
  }
  slots_=Index(coefficients.size());
  std::vector<float> planes(std::size_t(3*terms)*std::max<Index>(slots_,1),0);
  for(Index s=0;s<slots_;++s) for(int k=0;k<3*terms;++k) planes[std::size_t(k)*slots_+s]=coefficients[s][k];
  coefficients_=dev.upload<float>(planes);
  for(int level=0;level<=hops;++level) {
    std::vector<float> values(std::size_t(terms)*cells);
    for(Index c=0;c<cells;++c) for(int k=0;k<terms;++k) values[std::size_t(k)*cells+c]=float(geometry[level][c][k]);
    geometry_[level]=dev.upload<float>(values);
  }
  observed_host_=observed;observed_=dev.upload<Index>(observed);slot_=dev.upload<Index>(slot);mask_=dev.upload<Index>(mask);
  for(int k=0;k<2;++k) {moment_[k]=dev.alloc<float>(std::size_t(terms)*cells);lo_[k]=dev.alloc<float>(cells);hi_[k]=dev.alloc<float>(cells);}
}
bool WalkFaceReconstruction::observes(std::span<const BoundaryObservation> observations) const {
  if(observations.size()!=observed_host_.size()) return false;
  for(std::size_t b=0;b<observations.size();++b)
    if((observations[b]==BoundaryObservation::FaceValue)!=(observed_host_[b]!=0)) return false;
  return true;
}
std::shared_ptr<const WalkFaceReconstruction> walk_reconstruction(std::shared_ptr<const WalkFaceReconstruction> shared,
    Device& dev,const PolyMesh& mesh,const MeshGeometry& geometry,const DeviceMesh& device_mesh,const Stencil& stencil,
    std::span<const BoundaryObservation> observations) {
  if(shared && shared->observes(observations)) return shared;
  auto built=std::make_shared<const WalkFaceReconstruction>(dev,mesh,geometry,device_mesh,stencil,observations);
  return built->corrected_faces()>0?built:nullptr;
}
void WalkFaceReconstruction::moments(const Buffer<float>& field,const Buffer<float>& boundary) const {
  const auto& m=mesh_;
  kern::fv_moment_seed(dev_,m.n_cells,m.cell_face_offsets,stencil_.adjacent,observed_,m.owner_to_face.x,m.owner_to_face.y,m.owner_to_face.z,
    field,boundary,moment_[0],lo_[0],hi_[0],m.n_cells,m.n_internal_faces);
  for(int level=0;level<hops;++level) {
    const int from=level%2,to=1-from;
    kern::fv_moment_hop(dev_,m.n_cells,m.cell_face_offsets,m.cell_faces,stencil_.adjacent,m.cell_face_sign,
      m.owner_to_neighbour.x,m.owner_to_neighbour.y,m.owner_to_neighbour.z,field,geometry_[level],
      moment_[from],lo_[from],hi_[from],moment_[to],lo_[to],hi_[to],m.n_cells);
  }
  static_assert(hops%2==1,"moments() leaves the result in the second buffer");
}
void WalkFaceReconstruction::diffusion_correction(const Buffer<float>& field,const Buffer<float>& boundary,Buffer<float>& out) const {
  const auto& m=mesh_;moments(field,boundary);
  kern::fv_moment_diffusion(dev_,m.n_faces,m.owner,m.neighbour,slot_,m.owner_to_face.x,m.owner_to_face.y,m.owner_to_face.z,
    m.neighbour_to_face.x,m.neighbour_to_face.y,m.neighbour_to_face.z,field,geometry_[hops],moment_[1],coefficients_,out,
    m.n_cells,m.n_internal_faces,slots_);
}
void WalkFaceReconstruction::momentum_values(const Buffer<float>& field,const Buffer<float>& boundary,
    const Buffer<float>& diagonal,Buffer<float>& out,const Buffer<float>* cell_limiter) const {
  const auto& m=mesh_;moments(field,boundary);
  kern::fv_moment_momentum(dev_,m.n_internal_faces,m.owner,m.neighbour,slot_,m.owner_to_face.x,m.owner_to_face.y,m.owner_to_face.z,
    m.neighbour_to_face.x,m.neighbour_to_face.y,m.neighbour_to_face.z,field,geometry_[hops],moment_[1],coefficients_,
    lo_[1],hi_[1],diagonal,cell_limiter?*cell_limiter:diagonal,out,m.n_cells,m.n_internal_faces,slots_,cell_limiter?1:0);
}
} // namespace kutta::fv
