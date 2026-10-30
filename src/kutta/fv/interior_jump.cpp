#include "kutta/fv/interior_jump.hpp"
#include <algorithm>
#include <cmath>
#include <set>
#include <map>
#include "kutta/mesh/topology.hpp"
#include <limits>
#include "kernels.hpp"
#include "kutta/mesh/detail.hpp"

namespace kutta::fv {
using namespace detail;
InteriorJumps::InteriorJumps(Device& dev,const PolyMesh& mesh,const MeshGeometry& g,
    std::span<const InteriorCondition> conditions,double rho,double mu)
    :dev_(dev),n_(mesh.n_cells),nf_(mesh.n_faces()),rho_(float(rho)),mu_(float(mu)) {
  if(conditions.empty()) return;
  std::vector<OrientedFace> selected;
  for(const auto& condition:conditions) {
    const auto it=std::find_if(mesh.internal_face_zones.begin(),mesh.internal_face_zones.end(),[&](const auto& z){return z.source_zone_id==condition.source_zone_id;});
    KUTTA_REQUIRE(it!=mesh.internal_face_zones.end() && !it->faces.empty(),"jump: internal face zone absent from mesh");
    selected.insert(selected.end(),it->faces.begin(),it->faces.end());
  }
  active_faces_=index(selected.size());const auto stride=active_faces_;
  std::vector<Index> kind(stride,0),counts(stride,0),faces;std::vector<float> p(std::size_t(stride)*37,0),area;
  std::map<Index,Index> face_index;
  for(auto [f,direction]:selected) {
    KUTTA_REQUIRE(f>=0 && f<mesh.n_internal_faces() && (direction==1 || direction==-1) && face_index.emplace(f,index(faces.size())).second,
      "jump: invalid or overlapping oriented face");
    faces.push_back(f);area.push_back(float(norm(g.face_area[f])));
  }
  std::set<int> zones;
  for(const auto& c:conditions) {
    active_=true;KUTTA_REQUIRE(zones.insert(c.source_zone_id).second,"jump: duplicate face-zone condition");
    const auto it=std::find_if(mesh.internal_face_zones.begin(),mesh.internal_face_zones.end(),[&](const auto& z){return z.source_zone_id==c.source_zone_id;});
    KUTTA_REQUIRE(it!=mesh.internal_face_zones.end() && !it->faces.empty(),"jump: internal face zone absent from mesh");
    for(auto [f,direction]:it->faces) {
      const auto a=face_index.at(f);
      if(const auto* fan=std::get_if<FanCurve>(&c.condition)) {
        KUTTA_REQUIRE(!fan->coefficients.empty() && fan->coefficients.size()<=32,"fan: polynomial must have 1 to 32 coefficients");
        kind[a]=1;counts[a]=Index(fan->coefficients.size());p[a]=float(fan->reverse?-direction:direction);
        p[3*stride+a]=0;p[4*stride+a]=std::numeric_limits<float>::max();
        if(fan->velocity_range) {
          const auto r=*fan->velocity_range;
          KUTTA_REQUIRE(std::isfinite(r[0]) && std::isfinite(r[1]) && r[0]>=0 && r[1]>=r[0],"fan: invalid polynomial velocity range");
          p[3*stride+a]=float(r[0]);p[4*stride+a]=float(r[1]);
        }
        for(std::size_t j=0;j<fan->coefficients.size();++j) p[(5+j)*stride+a]=float(fan->coefficients[j]);
      } else {
        const auto& porous=std::get<PorousJump>(c.condition);
        KUTTA_REQUIRE(std::isfinite(porous.permeability) && porous.permeability>0 && std::isfinite(porous.thickness) && porous.thickness>=0 &&
          std::isfinite(porous.inertial) && porous.inertial>=0,"porous-jump: positive permeability and nonnegative thickness/C2 required");
        kind[a]=2;p[stride+a]=float(porous.thickness/porous.permeability);p[2*stride+a]=float(porous.thickness*porous.inertial);
      }
    }
  }
  for(auto v:p) KUTTA_REQUIRE(std::isfinite(v),"jump: coefficients overflow float storage");
  const auto cell_faces=build_cell_faces(mesh);
  std::map<Index,std::vector<std::pair<Index,Index>>> rows;
  for(auto f:faces) {rows.try_emplace(mesh.owner[f]);rows.try_emplace(mesh.neighbour[f]);}
  for(auto& [c,row]:rows) for(Index j=cell_faces.offsets[c];j<cell_faces.offsets[c+1];++j) {
    auto found=face_index.find(cell_faces.faces[j]);
    if(found!=face_index.end()) row.emplace_back(j,(cell_faces.signs[j]>0?1:-1)*(found->second+1));
  }
  std::vector<Index> cells,offsets{0},incidence,entries;
  for(const auto& [c,row]:rows) {
    cells.push_back(c);for(auto [j,a]:row) {incidence.push_back(j);entries.push_back(a);}
    offsets.push_back(index(entries.size()));
  }
  active_cells_=index(cells.size());
  faces_=dev.upload<Index>(faces);cells_=dev.upload<Index>(cells);offsets_=dev.upload<Index>(offsets);
  incidence_=dev.upload<Index>(incidence);entries_=dev.upload<Index>(entries);
  kinds_=dev.upload<Index>(kind);counts_=dev.upload<Index>(counts);parameters_=dev.upload<float>(p);area_=dev.upload<float>(area);
  jumps_=dev.alloc<float>(stride);slopes_=dev.alloc<float>(stride);previous_=dev.alloc<float>(stride);
  dev.fill(jumps_,0);dev.fill(slopes_,0);dev.fill(previous_,0);
}
void InteriorJumps::update(const Buffer<float>& mass) {
  if(active_) kern::fv_jump_update(dev_,active_faces_,faces_,kinds_,counts_,parameters_,area_,mass,jumps_,slopes_,previous_,active_faces_,rho_,mu_);
}
void InteriorJumps::gradient(const DeviceMesh&,const BoundaryLeastSquares& lsq,DeviceVec3& g) const {
  if(active_) kern::fv_jump_pressure_gradient(dev_,active_cells_,cells_,offsets_,incidence_,entries_,
    lsq.rows().x,lsq.rows().y,lsq.rows().z,lsq.inverse,jumps_,g.x,g.y,g.z,n_);
}
void InteriorJumps::momentum(const DeviceMesh& m,DeviceVec3& source) const {
  if(active_) kern::fv_jump_momentum(dev_,active_cells_,cells_,offsets_,entries_,faces_,m.area.x,m.area.y,m.area.z,m.volume,
    jumps_,source.x,source.y,source.z);
}
void InteriorJumps::flux(const DeviceMesh& m,const Stencil& s,const Buffer<float>& diagonal,
    Buffer<float>& mass,Buffer<float>& coefficient) const {
  if(active_) kern::fv_jump_flux(dev_,active_faces_,faces_,m.owner,m.neighbour,m.volume,diagonal,s.diffusion,
    jumps_,slopes_,previous_,mass,coefficient,rho_);
}
} // namespace kutta::fv
