#include "kutta/fv/face_projection.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <set>
#include <limits>
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
