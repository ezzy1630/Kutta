#include "kutta/setup/fluent_physics.hpp"
#include <algorithm>
#include <cmath>
#include <functional>
#include <map>
#include <limits>
#include "kutta/mesh/detail.hpp"

namespace kutta::setup {
using namespace detail;
namespace {
double scalar(const std::optional<fluent::ScalarValue>& v,const std::string& key) {
  const auto n=v?(v->constant?v->constant:v->resolved):std::nullopt;
  KUTTA_REQUIRE(n && std::isfinite(*n),"Fluent physics: "+key+" requires a finite constant or resolved expression");return *n;
}
Vec3d vector(const fluent::ScalarVector& v,const std::string& key) {
  return {scalar(v[0],key+"/x"),scalar(v[1],key+"/y"),scalar(v[2],key+"/z")};
}
bool flag(const fluent::RawSettings& raw,const std::string& key) {
  const auto* v=raw.value(key);
  KUTTA_REQUIRE(v && v->kind==fluent::SExpr::Kind::boolean,"Fluent physics: "+key+" requires an explicit boolean");return v->boolean;
}
void neutral_scalar(const fluent::RawSettings& raw,const std::string& key,double expected) {
  const auto* v=raw.value(key);
  KUTTA_REQUIRE(v && v->number()==expected,"Fluent physics: "+key+" requires "+std::to_string(expected));
}
fv::Rotation rotation(const fluent::ScalarVector& origin,const fluent::ScalarVector& direction,
    const std::optional<fluent::ScalarValue>& speed,const std::string& key) {
  const auto axis=vector(direction,key+"/axis");
  KUTTA_REQUIRE(norm(axis)>0,key+": nonzero rotation axis required");
  return {vector(origin,key+"/origin"),axis*(scalar(speed,key+"/omega")/norm(axis))};
}
}
fv::Wall fluent_wall(const fluent::BoundaryZone& z) {
  const auto key="zone/"+std::to_string(z.id);fv::Wall wall;
  KUTTA_REQUIRE(z.wall_motion && z.relative_motion,key+": explicit wall motion and relative? required");
  wall.relative_to_zone=*z.relative_motion;
  if(z.wall_motion==fluent::WallMotion::rotating) wall.rotation=rotation(z.axis_origin,z.axis_direction,z.angular_velocity,key);
  if(z.wall_motion==fluent::WallMotion::translating) {
    if(z.velocity_spec==fluent::VelocitySpec::components) wall.velocity=vector(z.velocity,key+"/velocity");
    else {
      KUTTA_REQUIRE(z.velocity_spec==fluent::VelocitySpec::magnitude_direction,key+": unsupported wall velocity representation");
      const auto d=vector(z.direction,key+"/direction");KUTTA_REQUIRE(norm(d)>0,key+": nonzero wall direction required");
      wall.velocity=d*(scalar(z.velocity_magnitude,key+"/vmag")/norm(d));
    }
  }
  return wall;
}
std::vector<fv::CellZonePhysics> fluent_zone_physics(const fluent::FluentCase& c) {
  std::vector<fv::CellZonePhysics> result;
  std::map<std::int64_t,const fluent::CellZone*> zones;std::map<std::int64_t,fv::RigidFrame> frames;
  std::map<std::int64_t,int> state;
  for(const auto& z:c.cell_zones) zones.emplace(z.id,&z);
  const bool moving=std::any_of(c.cell_zones.begin(),c.cell_zones.end(),[](const auto& z){return z.frame_motion==true;});
  const bool porous=std::any_of(c.cell_zones.begin(),c.cell_zones.end(),[](const auto& z){return z.porous==true;});
  if(moving) KUTTA_REQUIRE(flag(c.variables,"solve-absolute-velocities?"),"MRF: only Fluent's absolute-velocity formulation is implemented");
  if(porous) KUTTA_REQUIRE(!flag(c.variables,"porous/physical-velocity-form?"),"porous: only Fluent's superficial-velocity formulation is implemented");
  std::function<fv::RigidFrame(std::int64_t)> resolve=[&](std::int64_t id) {
    KUTTA_REQUIRE(zones.contains(id),"MRF: relative frame refers to an absent cell zone");
    KUTTA_REQUIRE(state[id]!=1,"MRF: cycle in relative frame references");
    if(state[id]==2) return frames.at(id);
    state[id]=1;const auto& z=*zones.at(id);fv::RigidFrame f;
    if(z.frame_motion==true) {
      const auto key="zone/"+std::to_string(id);
      KUTTA_REQUIRE(z.motion_spec==0,key+": only constant frame motion (motion-spec=0) is supported");
      const auto* udf=z.raw.value("mrf-udf-zmotion-name");
      KUTTA_REQUIRE(!udf || udf->atom_text()=="none",key+": UDF frame motion is unsupported");
      f={vector(z.frame_velocity,key+"/frame-velocity"),rotation(z.axis_origin,z.axis_direction,z.angular_velocity,key)};
      KUTTA_REQUIRE(z.relative_to_zone,key+": mrf-relative-to-thread is required");
      if(*z.relative_to_zone!=-1) {
        const auto parent=resolve(*z.relative_to_zone);
        if(!flag(c.variables,"reference-frame/relative/origin-absolute?")) f.rotation.origin=f.rotation.origin+parent.rotation.origin;
        f=fv::compose_frames(parent,f);
      }
    }
    state[id]=2;frames.emplace(id,f);return f;
  };
  for(const auto& z:c.cell_zones) {
    KUTTA_REQUIRE(z.frame_motion && z.porous,"Fluent physics: explicit mrf-motion? and porous? required");
    if(!*z.frame_motion && !*z.porous) continue;
    fv::CellZonePhysics p;KUTTA_REQUIRE(z.id>=0 && z.id<=std::numeric_limits<int>::max(),"physics: zone ID outside mesh identity range");p.source_zone_id=int(z.id);
    if(*z.frame_motion) p.frame=resolve(z.id);
    if(*z.porous) {
      const auto key="zone/"+std::to_string(z.id);
      KUTTA_REQUIRE(z.relative_resistance,key+": rel-vel-resistance? is required");
      KUTTA_REQUIRE(!flag(z.raw,"alt-inertial-form?"),key+": alternative inertial formulation is unsupported");
      if(z.raw.entry("dir-spec-cond")) neutral_scalar(z.raw,"dir-spec-cond",0);
      neutral_scalar(z.raw,"c0",0); // power-law source
      neutral_scalar(z.raw,"c1",0);
      fv::PorousResistance r;
      r.viscous=vector(z.viscous_resistance,key+"/viscous-resistance");r.inertial=vector(z.inertial_resistance,key+"/inertial-resistance");
      r.direction1=vector(z.porous_directions[0],key+"/direction-1");r.direction2=vector(z.porous_directions[1],key+"/direction-2");
      r.porosity=scalar(z.porosity,key+"/porosity");r.relative_velocity=*z.relative_resistance;
      (void)fv::porous_basis(r);p.porous=r;
      // A non-unit relative viscosity changes diffusion, beyond the basic model.
      const auto* ratio=z.raw.entry("viscosity-ratio");
      if(ratio) {
        auto value=fluent::entry_value(*ratio);
        if(value && value->kind==fluent::SExpr::Kind::list) value=fluent::entry_value(*value);
        KUTTA_REQUIRE(value && value->number()==1,"porous: relative viscosity models are unsupported");
      }
    }
    result.push_back(p);
  }
  return result;
}
} // namespace kutta::setup
