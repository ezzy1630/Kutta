#include "kutta/fv/zone_physics.hpp"

#include <algorithm>
#include <cmath>
#include <set>
#include <map>

#include "kernels.hpp"
#include "kutta/mesh/detail.hpp"

namespace kutta::fv {
using namespace detail;
Vec3d RigidFrame::velocity(Vec3d p) const {
  return translation+cross(rotation.angular_velocity,p-rotation.origin);
}
RigidFrame compose_frames(const RigidFrame& parent,const RigidFrame& child) {
  return {parent.velocity(child.rotation.origin)+child.translation,
          {child.rotation.origin,parent.rotation.angular_velocity+child.rotation.angular_velocity}};
}
std::array<Vec3d,3> porous_basis(const PorousResistance& p) {
  KUTTA_REQUIRE(finite(p.direction1) && finite(p.direction2) && norm(p.direction1)>0,
                "porous: finite nonzero Direction-1 required");
  auto a=p.direction1/norm(p.direction1);
  auto b=p.direction2-a*dot(a,p.direction2);
  KUTTA_REQUIRE(norm(b)>1e-12*norm(p.direction2),"porous: Direction-2 must not be parallel to Direction-1");
  b=b/norm(b);return {a,b,cross(a,b)};
}
namespace {
std::vector<const CellZonePhysics*> assignments(const PolyMesh& mesh,std::span<const CellZonePhysics> zones) {
  std::vector<const CellZonePhysics*> result(mesh.n_cells,nullptr);std::set<int> ids;
  for(const auto& z:zones) {
    KUTTA_REQUIRE(z.source_zone_id>=0 && ids.insert(z.source_zone_id).second,"physics: invalid or duplicate cell-zone identity");
    const auto it=std::find_if(mesh.cell_zones.begin(),mesh.cell_zones.end(),[&](const auto& m){return m.source_zone_id==z.source_zone_id;});
    KUTTA_REQUIRE(it!=mesh.cell_zones.end() && !it->cells.empty(),"physics: cell zone not present in mesh");
    for(auto c:it->cells) {
      KUTTA_REQUIRE(c>=0 && c<mesh.n_cells && !result[c],"physics: overlapping or invalid cell-zone cells");
      result[c]=&z;
    }
  }
  return result;
}
void validate_frame(const RigidFrame& f) {
  KUTTA_REQUIRE(finite(f.translation) && finite(f.rotation.origin) && finite(f.rotation.angular_velocity),"MRF: nonfinite rigid motion");
}
}
std::vector<RigidFrame> cell_frames(const PolyMesh& mesh,std::span<const CellZonePhysics> zones) {
  const auto assigned=assignments(mesh,zones);std::vector<RigidFrame> frames(mesh.n_cells);
  for(Index c=0;c<mesh.n_cells;++c) if(assigned[c] && assigned[c]->frame) {
    frames[c]=*assigned[c]->frame;validate_frame(frames[c]);
  }
  return frames;
}
ZonePhysics::ZonePhysics(Device& dev,const PolyMesh& mesh,const MeshGeometry& g,
    std::span<const CellZonePhysics> zones,double rho,double mu)
    :dev_(dev),n_(mesh.n_cells),nf_(mesh.n_faces()),rho_(float(rho)),mu_(float(mu)) {
  if(zones.empty()) return;
  const auto assigned=assignments(mesh,zones);zones_=index(zones.size());
  const RigidFrame stationary;
  const auto frame=[&](Index c)->const RigidFrame& {return assigned[c] && assigned[c]->frame?*assigned[c]->frame:stationary;};
  std::map<int,Index> zone_ids;
  std::vector<float> data(std::size_t(zones_)*20,0),flux;
  for(Index c=0;c<zones_;++c) {
    const auto& z=zones[c];zone_ids.emplace(z.source_zone_id,c);
    if(z.frame) {
      active_=moving_=true;validate_frame(*z.frame);const auto w=z.frame->rotation.angular_velocity;
      data[12*zones_+c]=float(w.x);data[13*zones_+c]=float(w.y);data[14*zones_+c]=float(w.z);
      data[15*zones_+c]=float(z.frame->translation.x);data[16*zones_+c]=float(z.frame->translation.y);data[17*zones_+c]=float(z.frame->translation.z);
    }
    if(!z.porous) continue;
    active_=true;const auto& p=*z.porous;
    relative_motion_=relative_motion_ || (z.frame && p.relative_velocity);
    KUTTA_REQUIRE(finite(p.viscous) && finite(p.inertial) && p.viscous.x>=0 && p.viscous.y>=0 && p.viscous.z>=0 &&
      p.inertial.x>=0 && p.inertial.y>=0 && p.inertial.z>=0 && std::isfinite(p.porosity) && p.porosity>=0 && p.porosity<=1,
      "porous: nonnegative resistances and porosity in [0,1] required");
    const auto basis=porous_basis(p);
    const std::array<double,3> d{p.viscous.x,p.viscous.y,p.viscous.z},a{p.inertial.x,p.inertial.y,p.inertial.z};
    for(int t=0;t<2;++t) {
      std::array<double,6> tensor{};const auto& values=t?a:d;
      for(int j=0;j<3;++j) {const auto e=basis[j];const auto v=values[j];
        tensor[0]+=v*e.x*e.x;tensor[1]+=v*e.x*e.y;tensor[2]+=v*e.x*e.z;
        tensor[3]+=v*e.y*e.y;tensor[4]+=v*e.y*e.z;tensor[5]+=v*e.z*e.z;
      }
      for(int j=0;j<6;++j) data[(t*6+j)*zones_+c]=float(tensor[j]);
      data[(18+t)*zones_+c]=float(*std::max_element(values.begin(),values.end()));
    }
  }
  std::vector<Index> indices(n_,-1);
  std::vector<float> reference_velocity(relative_motion_?std::size_t(n_)*3:1,0);
  for(Index c=0;c<n_;++c) if(const auto* z=assigned[c]) {
    indices[c]=zone_ids.at(z->source_zone_id);
    if(relative_motion_ && z->porous && z->porous->relative_velocity) {
      const auto v=frame(c).velocity(g.cell_centroid[c]);reference_velocity[c]=float(v.x);
      reference_velocity[n_+c]=float(v.y);reference_velocity[2*n_+c]=float(v.z);
    }
  }
  if(moving_) flux.resize(nf_,0);
  for(Index f=0;moving_ && f<nf_;++f) {
    const auto p=mesh.owner[f];const auto vf=frame(p).velocity(g.face_centroid[f]);
    flux[f]=float(rho*dot(vf,g.face_area[f]));
    if(f<mesh.n_internal_faces()) {
      const auto vn=frame(mesh.neighbour[f]).velocity(g.face_centroid[f]);
      // A conformal MRF interface is stationary in both frames. Never average
      // incompatible normal mesh speeds into a fictitious mass flux.
      if(norm(vf-vn)>1e-12) {
        const auto scale=std::max({norm(vf),norm(vn),1.0})*norm(g.face_area[f]);
        KUTTA_REQUIRE(std::abs(dot(vf,g.face_area[f]))<1e-6*scale && std::abs(dot(vn,g.face_area[f]))<1e-6*scale,
          "MRF: zone interface must have zero normal frame velocity on both sides");
        flux[f]=0;
      }
    }
  }
  for(const auto& patch:mesh.patches) if(moving_ && patch.kind==PatchKind::Symmetry)
    for(Index f=patch.start;f<patch.start+patch.size;++f) {
      const auto v=frame(mesh.owner[f]).velocity(g.face_centroid[f]);
      KUTTA_REQUIRE(std::abs(dot(v,g.face_area[f]))<1e-6*std::max(norm(v),1.0)*norm(g.face_area[f]),
        "MRF: symmetry boundary must have zero normal frame velocity");
    }
  for(auto v:data) KUTTA_REQUIRE(std::isfinite(v),"physics: coefficients overflow float storage");
  for(auto v:flux) KUTTA_REQUIRE(std::isfinite(v),"MRF: frame flux overflows float storage");
  for(auto v:reference_velocity) KUTTA_REQUIRE(std::isfinite(v),"porous: reference velocity overflows float storage");
  zone_index_=dev.upload<Index>(indices);relative_velocity_=dev.upload<float>(reference_velocity);
  coefficients_=dev.upload<float>(data);if(moving_) frame_flux_=dev.upload<float>(flux);
}
void ZonePhysics::relative_flux(const Buffer<Index>& kind,Index internal,Buffer<float>& mass) const {
  if(moving_) kern::fv_frame_flux(dev_,nf_,frame_flux_,kind,mass,internal);
}
void ZonePhysics::momentum(int component,const DeviceVec3& u,const Buffer<float>& volume,
    Buffer<float>& diagonal,Buffer<float>& rhs) const {
  if(active_) kern::fv_zone_momentum(dev_,n_,zone_index_,coefficients_,relative_velocity_,u.x,u.y,u.z,volume,diagonal,rhs,n_,zones_,component,relative_motion_?1:0,rho_,mu_);
}
} // namespace kutta::fv
