#include "kutta/setup/kutta_case.hpp"
#include "kutta/setup/fluent.hpp"
#include "kutta/setup/settings.hpp"
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
