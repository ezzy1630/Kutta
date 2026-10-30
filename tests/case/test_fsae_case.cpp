#include <doctest/doctest.h>
#include <algorithm>
#include "kutta/setup/fluent_physics.hpp"
#include "kutta/fv/boundary_physics.hpp"
#include "kutta/mesh/detail.hpp"
#include "kutta/io/fluent_detail.hpp"
#include "kutta/mesh/topology.hpp"

using namespace kutta;
using namespace kutta::detail;
namespace f=fluent;
namespace {
f::FluentCase extracted(const std::string& threads,const std::string& variables="") {
  return f::extract_case({{"Rampant Variables",f::parse_sexprs("(37 ((case-config ((rp-3d? . #t) (rp-seg? . #t) (rp-unsteady? . #f) (rp-lam? . #t) (rp-kw? . #f) (rp-ke? . #f))) "+variables+"))")},
    {"Thread Variables",f::parse_sexprs(threads)}});
}
std::string frame(int id,int parent,double omega,const std::string& origin) {
  return "(39 ("+std::to_string(id)+" fluid rotor 1) ((mrf-motion? . #t) (porous? . #f) (motion-spec . 0) (mrf-relative-to-thread . "+std::to_string(parent)+") "
    "(mrf-omega . "+std::to_string(omega)+") (mrf-grid-x-vel . .25) (mrf-grid-y-vel . 0) (mrf-grid-z-vel . 0) "
    "(reference-frame-axis-origin-components "+origin+") (reference-frame-axis-direction-components 0 0 2)))";
}
}

TEST_CASE("FSAE Fluent wall representations retain SI radians and relative motion") {
  const auto c=extracted(R"(
    (39 (4 wall ground 1) ((motion-bc . 1) (rotating? . #f) (components? . #f) (relative? . #f)
      (vmag . 12) (ni . 3) (nj . 0) (nk . 0)))
    (39 (5 wall wheel 1) ((motion-bc . 1) (rotating? . #t) (relative? . #t) (omega . 2)
      (x-origin . 1) (y-origin . 0) (z-origin . 0) (ai . 0) (aj . 0) (ak . 3)))
  )");
  const auto ground=setup::fluent_wall(c.boundaries[0]),wheel=setup::fluent_wall(c.boundaries[1]);
  CHECK(ground.velocity.x==12);CHECK_FALSE(ground.relative_to_zone);
  CHECK(wheel.rotation->angular_velocity.z==2);CHECK(wheel.relative_to_zone);
  const auto velocity=fv::wall_velocity(wheel,{2,0,0},fv::RigidFrame{{1,0,0},{{},{0,0,3}}});
  CHECK(velocity.x==1);CHECK(velocity.y==8);
}

TEST_CASE("FSAE Fluent nested MRF frames resolve origins and reject cycles") {
  const std::string variables="(solve-absolute-velocities? . #t) (reference-frame/relative/origin-absolute? . #t)";
  auto c=extracted(frame(7,-1,2,"1 0 0")+frame(8,7,3,"0 2 0"),variables);
  const auto zones=setup::fluent_zone_physics(c);REQUIRE(zones.size()==2);
  const auto& child=*zones[1].frame;
  CHECK(child.rotation.angular_velocity.z==5);
  const auto expected=fv::RigidFrame{{.25,0,0},{{1,0,0},{0,0,2}}}.velocity({4,5,0})+
    fv::RigidFrame{{.25,0,0},{{0,2,0},{0,0,3}}}.velocity({4,5,0});
  CHECK(detail::norm(child.velocity({4,5,0})-expected)<1e-12);
  c=extracted(frame(7,8,2,"1 0 0")+frame(8,7,3,"0 2 0"),variables);
  CHECK_THROWS_WITH(setup::fluent_zone_physics(c),doctest::Contains("cycle"));
  c=extracted(frame(7,-1,2,"1 0 0"),"(solve-absolute-velocities? . #f)");
  CHECK_THROWS_WITH(setup::fluent_zone_physics(c),doctest::Contains("absolute-velocity"));
}

TEST_CASE("FSAE Fluent porous coordinates use orthogonalized directions and superficial velocity") {
  const std::string zone=R"((39 (7 fluid radiator 1) ((mrf-motion? . #f) (porous? . #t) (rel-vel-resistance? . #t)
    (porous-r/1 . 2) (porous-r/2 . 7) (porous-r/3 . 4) (porous-c/1 . 3) (porous-c/2 . 1) (porous-c/3 . 5)
    (direction-1/x . 2) (direction-1/y . 1) (direction-1/z . 0)
    (direction-2/x . 1) (direction-2/y . 3) (direction-2/z . 1)
    (porosity . .35) (alt-inertial-form? . #f) (c0 . 0) (c1 . 0) (viscosity-ratio (constant . 1)))))";
  auto c=extracted(zone,"(porous/physical-velocity-form? . #f)");const auto result=setup::fluent_zone_physics(c);
  REQUIRE(result.size()==1);CHECK(result[0].porous->porosity==.35);
  const auto basis=fv::porous_basis(*result[0].porous);CHECK(std::abs(detail::dot(basis[0],basis[1]))<1e-12);
  CHECK(detail::dot(basis[2],detail::cross(basis[0],basis[1]))==doctest::Approx(1));
  c=extracted(zone,"(porous/physical-velocity-form? . #t)");
  CHECK_THROWS_WITH(setup::fluent_zone_physics(c),doctest::Contains("superficial"));
}

TEST_CASE("FSAE Fluent pressure mass and porous jump inputs remain distinct") {
  auto c=extracted(R"(
    (39 (4 pressure-inlet inlet 1) ((p0 . 12) (p . 0) (frame-of-reference . 0) (direction-spec . 1) (flow-spec . 0) (coordinate-system . 0)))
    (39 (5 mass-flow-inlet rate 1) ((mass-flow . .2) (flow-spec . 0) (p . 0) (frame-of-reference . 0) (direction-spec . 1)))
    (39 (6 mass-flow-inlet flux 1) ((mass-flux . 2) (flow-spec . 1) (p . 0) (frame-of-reference . 0) (direction-spec . 1)))
    (39 (7 porous-jump radiator 1) ((alpha . .01) (dm . .2) (c2 . 4)))
    (39 (8 pressure-far-field farfield 1) ())
  )");
  REQUIRE(c.boundaries.size()==5);
  CHECK(c.boundaries[0].type==f::ZoneType::pressure_inlet);CHECK(c.boundaries[0].total_pressure->constant==12);
  CHECK(c.boundaries[1].mass_flow_spec==0);CHECK(c.boundaries[1].mass_flow->constant==.2);
  CHECK(c.boundaries[2].mass_flow_spec==1);CHECK(c.boundaries[2].mass_flux->constant==2);
  CHECK(c.boundaries[3].jump_permeability->constant==.01);CHECK(c.boundaries[3].jump_thickness->constant==.2);
  CHECK(c.boundaries[3].jump_inertial->constant==4);
  CHECK(std::any_of(c.issues.begin(),c.issues.end(),[](const auto& i){return i.path=="zone/8" && i.reason.find("compressible")!=std::string::npos;}));
  for(const auto& issue:c.issues) CHECK(issue.path!="zone/6/flow-spec");
}

TEST_CASE("FSAE Fluent import retains special interior identity and input normal through RCM") {
  const auto box=make_box_hex(2,1,1,{0,0,0},{1,1,.1});const auto geometry=compute_geometry(box);
  FluentRaw raw;raw.zones.emplace(7,FluentZone{"fluid","fluid"});raw.zones.emplace(71,FluentZone{"fan","fan"});
  for(Index p=0;p<box.n_points();++p) raw.nodes.emplace(p+1,box.points[p]);
  for(Index c=0;c<box.n_cells;++c) raw.cells.emplace(c+1,7);
  for(Index f=0;f<box.n_faces();++f) {
    const auto face=nodes(box,f);std::vector<Index> source_nodes;
    for(auto p:face) source_nodes.push_back(p+1);std::reverse(source_nodes.begin(),source_nodes.end());
    const auto zone=f<box.n_internal_faces()?71:100;
    raw.faces.emplace(f+1,FluentFace{source_nodes,box.owner[f]+1,f<box.n_internal_faces()?box.neighbour[f]+1:0,zone});
  }
  raw.zones.emplace(100,FluentZone{"walls","wall"});
  auto imported=assemble_fluent(raw,1).mesh;
  REQUIRE(imported.internal_face_zones.size()==1);CHECK(imported.internal_face_zones[0].source_zone_id==71);
  renumber_cells_rcm(imported);const auto g=compute_geometry(imported);Vec3d normal{};
  for(auto [face,direction]:imported.internal_face_zones[0].faces) normal=normal+g.face_area[face]*double(direction);
  CHECK(norm(normal+geometry.face_area[0])<1e-12);
}
