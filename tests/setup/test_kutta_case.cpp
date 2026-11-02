#include "kutta/setup/kutta_case.hpp"
#include "kutta/setup/fluent.hpp"
#include "kutta/setup/settings.hpp"
#include "kutta/setup/fluent_physics.hpp"
#include "kutta/io/mesh_io.hpp"
#include "../io/test_support.hpp"
#include <doctest/doctest.h>
#include <algorithm>

namespace {
using namespace kutta;
using namespace kutta::setup;
const std::filesystem::path box_mesh = KUTTA_SOURCE_DIR "/tests/data/meshes/fsae_box.msh";

// The issues a rejected file reports, as "path: reason" lines.
std::vector<std::string> rejection(std::string_view text) {
  try {
    (void)parse_kutta_case(text);
  } catch (const RejectedSettings& e) {
    std::vector<std::string> out;
    for (const auto& i : e.issues()) out.push_back(i.path + ": " + i.reason);
    return out;
  }
  return {};
}
bool mentions(const std::vector<std::string>& lines, std::string_view text) {
  return std::any_of(lines.begin(), lines.end(), [&](const std::string& l) { return l.find(text) != std::string::npos; });
}
}  // namespace

TEST_CASE("Kutta case files reject unknown properties, wrong types and missing requirements with JSON paths") {
  const auto minimal = R"({"kutta_case": 1, "mesh": "box.msh", "models": {"viscous": "laminar"}, "boundaries": {"shell": {"type": "wall"}}})";
  CHECK(rejection(minimal).empty());
  const auto c = parse_kutta_case(minimal);
  CHECK(c.models.viscous == ViscousModel::laminar);
  CHECK(c.mesh == "box.msh");
  CHECK(std::holds_alternative<WallInput>(c.boundaries.at("shell")));

  auto issues = rejection(R"({"kutta_case": 1, "mesh": "box.msh", "models": {"viscous": "laminar", "viscus": 1},
    "boundaries": {"inlet": {"type": "velocity-inlet", "sped": 3}}})");
  CHECK(mentions(issues, "models/viscus: unknown property (did you mean \"viscous\"?)"));
  CHECK(mentions(issues, "boundaries/inlet/sped: unknown property (did you mean \"speed\"?)"));

  issues = rejection(R"({"kutta_case": 1, "mesh": "box.msh", "models": {"viscous": "k-omega"},
    "material": {"density": "1.2"}, "reference_values": {"area": -1}})");
  CHECK(mentions(issues, "models/viscous: \"k-omega\" is not one of: laminar, k-omega-sst, k-epsilon-standard"));
  CHECK(mentions(issues, "material/density: expected a number"));
  CHECK(mentions(issues, "reference_values/area: must be greater than 0"));

  issues = rejection(R"({"mesh": ""})");
  CHECK(mentions(issues, "kutta_case: required"));
  CHECK(mentions(issues, "mesh: required"));
  CHECK(mentions(issues, "models: required"));
  CHECK(mentions(issues, "boundaries: required"));
  CHECK(mentions(rejection(R"({"kutta_case": 2, "mesh": "m.msh", "models": {"viscous": "laminar"}})"), "unsupported version"));
  CHECK(mentions(rejection("{not json"), "invalid JSON"));
  // Settings that only make sense for another model are errors, not ignored.
  CHECK(mentions(rejection(R"({"kutta_case": 1, "mesh": "m.msh", "models": {"viscous": "laminar", "transition": "intermittency"}})"),
                 "models/transition"));
  CHECK(mentions(rejection(R"({"kutta_case": 1, "mesh": "m.msh", "models": {"viscous": "laminar", "k_omega_sst": {}}})"),
                 "applies only to viscous k-omega-sst"));
}

TEST_CASE("The FSAE starter suggests roles from zone names and lowers to a case Kutta can run") {
  const auto mesh = read_mesh(box_mesh);
  const auto info = mesh_info(mesh);
  CHECK(info.dimension == 3);
  CHECK_FALSE(info.polyhedral);
  std::map<std::string, std::string> roles;
  for (const auto& z : info.zones)
    if (z.kind == MeshZoneInfo::Kind::boundary) roles[z.name] = std::string(role_name(suggest_role(z).role));
  // The mesher exported every zone as a wall; names decide the roles.
  CHECK(roles == std::map<std::string, std::string>{{"front_wing", "wall"}, {"ground", "moving-ground"}, {"inlet", "velocity-inlet"},
                                                    {"outlet", "pressure-outlet"}, {"side_left", "symmetry"}, {"side_right", "symmetry"},
                                                    {"top", "symmetry"}});
  const auto wing = std::find_if(info.zones.begin(), info.zones.end(), [](const auto& z) { return z.name == "front_wing"; });
  CHECK(suggest_role(*wing).group == "Front wing");

  const auto starter = fsae_starter(info, "fsae_box.msh", {.speed = 20, .roles = {}});
  // The text form round-trips exactly: parse(text(c)) == c and text(parse(text)) == text.
  const auto text = kutta_case_text(starter);
  CHECK(parse_kutta_case(text) == starter);
  CHECK(kutta_case_text(parse_kutta_case(text)) == text);
  CHECK(text.find("\"speed\": 20") != std::string::npos);

  const auto lowered = lower_kutta_case(starter, info);
  auto retagged = mesh;
  apply_boundary_kinds(starter, retagged);
  const auto setup = fluent_setup(lowered, retagged);
  CHECK(setup.approximations.empty());
  const auto forces = force_definitions(lowered, retagged);
  REQUIRE(forces.size() == 2);
  for (const auto& f : forces) CHECK(f.patches == std::vector<std::string>{"front_wing"});

  // Overriding a role moves the zone, and the settings inventory sees it.
  auto as_wall = fsae_starter(info, "fsae_box.msh", {.speed = 20, .roles = {{"ground", ZoneRole::wall}}});
  CHECK(std::get<WallInput>(as_wall.boundaries.at("ground")).motion != WallInput::Motion::moving);
  RunControls run;
  auto lowered_wall = lower_kutta_case(as_wall, info);
  const auto inventory = settings_inventory(lowered_wall, run);
  CHECK(std::any_of(inventory.begin(), inventory.end(), [](const auto& s) { return s.path == "boundary/inlet/velocity-magnitude"; }));
}

TEST_CASE("Lowering checks a Kutta case against its mesh") {
  const auto info = mesh_info(read_mesh(box_mesh));
  auto c = fsae_starter(info, "fsae_box.msh", {});
  c.boundaries.erase("top");
  c.boundaries["tpo"] = SymmetryInput{};
  try {
    (void)lower_kutta_case(c, info);
    FAIL("expected RejectedSettings");
  } catch (const RejectedSettings& e) {
    std::vector<std::string> lines;
    for (const auto& i : e.issues()) lines.push_back(i.path + ": " + i.reason);
    CHECK(mentions(lines, "boundaries/tpo: the mesh has no zone named \"tpo\" (did you mean \"top\"?)"));
    CHECK(mentions(lines, "boundaries/top: required: mesh boundary zone \"top\""));
  }
}

TEST_CASE("A Kutta case file finds its mesh relative to itself and its identity covers both files") {
  MeshTemp t;
  std::filesystem::create_directory(t.path("meshes"));
  std::filesystem::copy_file(box_mesh, t.path("meshes/box.msh"));
  const auto info = mesh_info(read_mesh(box_mesh));
  mesh_text(t.path("case.kutta.json"), kutta_case_text(fsae_starter(info, "meshes/box.msh", {})));
  CHECK(is_kutta_case(t.path("case.kutta.json")));
  CHECK_FALSE(is_kutta_case(t.path("meshes/box.msh")));
  const auto loaded = load_kutta_case(t.path("case.kutta.json"));
  CHECK(loaded.mesh_path == (t.directory / "meshes/box.msh").lexically_normal());
  CHECK(loaded.mesh.n_cells == 70);
  CHECK(kutta_case_identity("a", "b") != kutta_case_identity("a", "c"));
  CHECK(kutta_case_identity("a", "b") != kutta_case_identity("b", "a"));
}

TEST_CASE("Kutta physics inputs round trip and lower to the Fluent boundary and cell-zone fields") {
  auto mesh = read_mesh(box_mesh);
  mesh.internal_face_zones.push_back({"screen", "porous-jump", 3, {{0,1}}});
  const auto info = mesh_info(mesh);
  auto c = fsae_starter(info, "fsae_box.msh", {});
  const auto inputs = parse_kutta_case(R"({"kutta_case":1,"mesh":"fsae_box.msh",
    "models":{"viscous":"k-omega-sst","transition":"intermittency",
      "k_omega_sst":{"cross_diffusion_clip":123},"intermittency":{"ctu1":90,"cpg2":-0.7,"sigma":1.2}},
    "cell_zones":{"fluid":{"frame_motion":{"velocity":[1,2,3],"angular_speed":4,"axis_origin":[0.1,0.2,0.3],"axis_direction":[0,2,0]},
      "porous":{"viscous_resistance":[30,200,400],"inertial_resistance":[4,20,40],"porosity":0.35,"relative_velocity":false,
      "direction_1":[1,1,0],"direction_2":[-1,1,0]}}},
    "pressure_reference_point":[0.2,0.1,0],
    "boundaries":{"inlet":{"type":"pressure-inlet","total_pressure":150,"direction":[2,0,0],"gauge_pressure":3,
      "turbulence":{"intensity":0.02},"intermittency":0.4},
      "front_wing":{"type":"wall","motion":"rotating","angular_speed":-60,"axis_origin":[0.2,0,0.2],"axis_direction":[0,1,0],"relative":false},
      "screen":{"type":"porous-jump","permeability":0.001,"thickness":0.02,"inertial_resistance":4}}})");
  c.models = inputs.models;
  c.cell_zones = inputs.cell_zones;
  c.pressure_reference_point = inputs.pressure_reference_point;
  for (const auto& [name, b] : inputs.boundaries) c.boundaries[name] = b;
  CHECK(parse_kutta_case(kutta_case_text(c)) == c);
  const auto lowered = lower_kutta_case(c, info);
  apply_boundary_kinds(c, mesh);
  const auto problem = fluent_setup(lowered, mesh);
  REQUIRE(problem.settings.cell_zone_physics.size() == 1);
  const auto& cell = problem.settings.cell_zone_physics.front();
  REQUIRE(cell.frame);
  CHECK(cell.frame->translation.x == 1);
  CHECK(cell.frame->rotation.angular_velocity.y == 4);
  CHECK(cell.frame->rotation.origin.z == 0.3);
  REQUIRE(cell.porous);
  CHECK(cell.porous->viscous.z == 400);
  CHECK(cell.porous->porosity == 0.35);
  CHECK_FALSE(cell.porous->relative_velocity);
  const auto inlet = std::find_if(problem.boundaries.begin(), problem.boundaries.end(), [](const auto& b) { return b.patch == "inlet"; });
  REQUIRE(inlet != problem.boundaries.end());
  CHECK(std::get<fv::PressureInlet>(inlet->condition).total_pressure == 150);
  REQUIRE(problem.settings.interior_conditions.size() == 1);
  CHECK(std::get<fv::PorousJump>(problem.settings.interior_conditions[0].condition).permeability == 0.001);
  CHECK(lowered.viscous.transition_coefficients.at("trans-one-eqn-ctu1") == 90);
  CHECK(lowered.viscous.transition_coefficients.at("trans-one-eqn-cpg2") == -0.7);
  CHECK(lowered.viscous.transition_coefficients.at("trans-one-eqn-ctu2") == 1000);
  CHECK(lowered.viscous.transition_coefficients.at("trans-intermit/min-limit") == 1e-7);
  CHECK(lowered.variables.value("sst-limit-cross-diff-to-n-times-w-destruct")->number() == 123);
  auto converted = kutta_case_from_fluent(lowered, info, c.mesh);
  CHECK(converted.omitted.empty());
  const auto again = lower_kutta_case(parse_kutta_case(kutta_case_text(converted.kutta)), info);
  CHECK(again.boundaries.size() == lowered.boundaries.size());
  CHECK(again.cell_zones[0].frame_velocity == lowered.cell_zones[0].frame_velocity);
  CHECK(again.viscous.transition_coefficients == lowered.viscous.transition_coefficients);
}

TEST_CASE("Kutta pressure mass-flow porous and frame defaults match Fluent defaults") {
  auto mesh = read_mesh(box_mesh);
  mesh.internal_face_zones.push_back({"screen", "porous-jump", 3, {{0,1}}});
  const auto info = mesh_info(mesh);
  auto c = fsae_starter(info, "fsae_box.msh", {});
  c.boundaries["inlet"] = PressureInletInput{};
  c.cell_zones["fluid"] = CellZoneInput{{}, FrameInput{}, PorousInput{}};
  c.boundaries["screen"] = PorousJumpInput{};
  auto lowered = lower_kutta_case(c, info);
  const auto inlet = std::find_if(lowered.boundaries.begin(), lowered.boundaries.end(), [](const auto& z) { return z.name == "inlet"; });
  REQUIRE(inlet != lowered.boundaries.end());
  CHECK(inlet->total_pressure->constant == 0);
  CHECK(inlet->direction_spec == 1);
  CHECK(inlet->frame_of_reference == 0);
  CHECK(inlet->turbulence.at("turb-intensity").constant == 0.05);
  CHECK(inlet->turbulence.at("turb-viscosity-ratio").constant == 10);
  REQUIRE(lowered.cell_zones.size() == 1);
  CHECK(lowered.cell_zones[0].relative_to_zone == -1);
  CHECK(lowered.cell_zones[0].porosity->constant == 1);
  CHECK(lowered.cell_zones[0].viscous_resistance[0]->constant == 0);
  CHECK(lowered.cell_zones[0].axis_direction[2]->constant == 1);
  CHECK_FALSE(lowered.variables.value("reference-frame/relative/origin-absolute?")->boolean);
  const auto jump = std::find_if(lowered.boundaries.begin(), lowered.boundaries.end(), [](const auto& z) { return z.name == "screen"; });
  REQUIRE(jump != lowered.boundaries.end());
  CHECK(jump->jump_permeability->constant == 1e10);
  CHECK(jump->jump_thickness->constant == 0);
  CHECK(jump->jump_inertial->constant == 0);
  c.boundaries["inlet"] = MassFlowInletInput{};
  lowered = lower_kutta_case(c, info);
  const auto mass = std::find_if(lowered.boundaries.begin(), lowered.boundaries.end(), [](const auto& z) { return z.name == "inlet"; });
  REQUIRE(mass != lowered.boundaries.end());
  CHECK(mass->mass_flow_spec == 0);
  CHECK(mass->mass_flow->constant == 1);
  c.boundaries["inlet"] = MassFlowInletInput{MassFlowInletInput::Flux{0.2}, {}, {}, {}, {}};
  lowered = lower_kutta_case(c, info);
  const auto flux = std::find_if(lowered.boundaries.begin(), lowered.boundaries.end(), [](const auto& z) { return z.name == "inlet"; });
  REQUIRE(flux != lowered.boundaries.end());
  CHECK(flux->mass_flow_spec == 1);
  CHECK(flux->mass_flux->constant == 0.2);
}

TEST_CASE("New Kutta physics rejects wrong types unknown keys and invalid combinations at JSON paths") {
  const auto issues = rejection(R"({"kutta_case":1,"mesh":"m.msh","models":{"viscous":"laminar"},
    "boundaries":{"inlet":{"type":"pressure-inlet","total_presure":3,"direction":[0,0,0]},
      "mass":{"type":"mass-flow-inlet","mass_flow_rate":1,"mass_flux":2},
      "jump":{"type":"porous-jump","permeability":0,"thickness":"1"}},
    "cell_zones":{"fluid":{"frame_motion":{"angular_speed":"fast","axis_direcion":[0,0,1]},
      "porous":{"porosity":2,"viscous_resistance":[1,-2,3],"direction_1":[1,0,0],"direction_2":[2,0,0]}}}})");
  CHECK(mentions(issues, "boundaries/inlet/total_presure: unknown property (did you mean \"total_pressure\"?)"));
  CHECK(mentions(issues, "boundaries/inlet/direction: must be a nonzero vector"));
  CHECK(mentions(issues, "boundaries/mass/mass_flux: give either mass_flow_rate or mass_flux"));
  CHECK(mentions(issues, "boundaries/jump/permeability: must be greater than 0"));
  CHECK(mentions(issues, "boundaries/jump/thickness: expected a number"));
  CHECK(mentions(issues, "cell_zones/fluid/frame_motion/angular_speed: expected a number"));
  CHECK(mentions(issues, "axis_direcion: unknown property (did you mean \"axis_direction\"?)"));
  CHECK(mentions(issues, "cell_zones/fluid/porous/porosity: must be at most 1"));
  CHECK(mentions(issues, "cell_zones/fluid/porous/viscous_resistance: components must be at least 0"));
  CHECK(mentions(issues, "cell_zones/fluid/porous/direction_2: must not be parallel"));
  const auto model = rejection(R"({"kutta_case":1,"mesh":"m.msh","models":{"viscous":"k-omega-sst","transition":"intermittency",
    "k_omega_sst":{"cross_diffusion_clip":0},"intermittency":{"sigma":0,"minimum":1,"maximum":0.5,"csep":-1}},"boundaries":{}})");
  CHECK(mentions(model, "models/k_omega_sst/cross_diffusion_clip: must be greater than 0"));
  CHECK(mentions(model, "models/intermittency/sigma: must be greater than 0"));
  CHECK(mentions(model, "models/intermittency/minimum: requires 0 < minimum < maximum <= 1"));
  CHECK(mentions(model, "models/intermittency/csep: must be at least 0"));
}

TEST_CASE("Frame references use cell-zone names reject cycles and compose relative origins") {
  const MeshInfo info{3, false, {{1,"parent","fluid",MeshZoneInfo::Kind::cell,1},{2,"child","fluid",MeshZoneInfo::Kind::cell,1}}, 2,0,0};
  auto c = parse_kutta_case(R"({"kutta_case":1,"mesh":"m.msh","models":{"viscous":"laminar"},"boundaries":{},
    "relative_frame_origins_absolute":false,"cell_zones":{
      "parent":{"frame_motion":{"angular_speed":1,"axis_origin":[1,0,0]}},
      "child":{"frame_motion":{"relative_to":"parent","axis_origin":[0,2,0],"angular_speed":2}}}})");
  const auto f = fluent_zone_physics(lower_kutta_case(c, info));
  REQUIRE(f.size() == 2);
  CHECK(f[1].frame->rotation.origin.x == 1);
  CHECK(f[1].frame->rotation.origin.y == 2);
  CHECK(f[1].frame->rotation.angular_velocity.z == 3);
  CHECK(f[1].frame->translation.x == -2);
  c.cell_zones["parent"].frame_motion->relative_to = "child";
  CHECK_THROWS_WITH(lower_kutta_case(c, info), doctest::Contains("cell_zones/child/frame_motion/relative_to: cycle"));
  c.cell_zones["parent"].frame_motion->relative_to.reset();
  c.cell_zones["child"].frame_motion->relative_to = "absent";
  CHECK_THROWS_WITH(lower_kutta_case(c, info), doctest::Contains("must name a mesh cell zone"));
}

TEST_CASE("Starter wheel geometry sets absolute rolling speed and preserves translating ground") {
  const auto info = mesh_info(read_mesh(box_mesh));
  StarterOptions o;
  o.speed = 20;
  o.wheels.emplace("front_wing", StarterOptions::Wheel{0.25, {0.2,0,0.25}, {0,-1,0}});
  const auto c = fsae_starter(info, "fsae_box.msh", o);
  const auto& w = std::get<WallInput>(c.boundaries.at("front_wing"));
  CHECK(w.motion == WallInput::Motion::rotating);
  CHECK(w.angular_speed == 80);
  CHECK(w.relative == false);
  CHECK(std::get<WallInput>(c.boundaries.at("ground")).speed == 20);
  o.wheels.at("front_wing").radius = 0;
  CHECK_THROWS_WITH(fsae_starter(info, "fsae_box.msh", o), doctest::Contains("wheel/front_wing/radius"));
}

TEST_CASE("Named fluid materials retain cell assignments and require uniform solver properties") {
  auto mesh = read_mesh(box_mesh);
  mesh.cell_zones.push_back({"second", "fluid", 20, {mesh.cell_zones[0].cells.back()}});
  mesh.cell_zones[0].cells.pop_back();
  const auto info = mesh_info(mesh);
  auto c = fsae_starter(info, "fsae_box.msh", {});
  c.materials.emplace("air-copy", KuttaCase::Material{{}, 1.225, 1.7894e-5});
  c.cell_zones["second"].material = "air-copy";
  CHECK(parse_kutta_case(kutta_case_text(c)) == c);
  const auto lowered = lower_kutta_case(c, info);
  REQUIRE(lowered.materials.size() == 2);
  CHECK(lowered.cell_zones[1].material == "air-copy");
  apply_boundary_kinds(c, mesh);
  CHECK_NOTHROW(fluent_setup(lowered, mesh));
  const auto converted = kutta_case_from_fluent(lowered, info, c.mesh);
  CHECK(converted.omitted.empty());
  CHECK(converted.kutta.cell_zones.at("second").material == "air-copy");
  CHECK(converted.kutta.materials.at("air-copy").density == 1.225);
  c.materials.at("air-copy").density = 2;
  CHECK_THROWS_WITH(fluent_setup(lower_kutta_case(c, info), mesh), doctest::Contains("constant uniform density/viscosity"));
  CHECK(mentions(rejection(R"({"kutta_case":1,"mesh":"m.msh","models":{"viscous":"laminar"},"boundaries":{},
    "materials":{"oil":{"density":1000,"viscocity":0.1}},"cell_zones":{"fluid":{"porous":{}}}})"),
    "materials/oil/viscocity: unknown property (did you mean \"viscosity\"?)"));
  CHECK(mentions(rejection(R"({"kutta_case":1,"mesh":"m.msh","models":{"viscous":"laminar"},"boundaries":{},
    "cell_zones":{"fluid":{"porous":{}}}})"), "cell_zones/fluid/porous/viscous_resistance: required"));
}
