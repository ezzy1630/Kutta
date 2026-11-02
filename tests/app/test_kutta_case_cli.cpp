// Kutta cases (mesh + *.kutta.json) through the real CLI: a case converted from
// a Fluent case solves to the same bytes as the Fluent case, and a starter
// written from a bare mesh runs.
#include "../backend_helpers.hpp"
#include "../fluent_examples.hpp"
#include "../generated_fixtures.hpp"
#include "../io/test_support.hpp"
#include "kutta/core/json_value.hpp"
#include "kutta/io/mesh_io.hpp"
#include <regex>
#include <set>
#include "kutta/setup/kutta_case.hpp"
#include "kutta/setup/fluent.hpp"
#include <iostream>
#include <doctest/doctest.h>
#include <cstdlib>
#include <string>
#include <vector>

namespace {
using namespace kutta;
std::string quote(const std::filesystem::path& path) {
#ifdef _WIN32
  return '"' + std::filesystem::path(path).make_preferred().string() + '"';
#else
  std::string out = "'";
  for (char c : path.string()) out += c == '\'' ? "'\\''" : std::string(1, c);
  return out + "'";
#endif
}
int kutta_cli(const std::string& arguments, const std::filesystem::path& log) {
  auto command = quote(KUTTA_CLI_PATH) + " " + arguments + " > " + quote(log) + " 2>&1";
#ifdef _WIN32
  command = '"' + command + '"';
#endif
  return std::system(command.c_str());
}
std::vector<std::string> lines(const std::filesystem::path& p) {
  std::vector<std::string> out;
  std::istringstream in(mesh_text(p));
  for (std::string l; std::getline(in, l);) out.push_back(l);
  return out;
}
// residuals.csv without its wall-clock column.
std::vector<std::string> residual_rows(const std::filesystem::path& p) {
  auto out = lines(p);
  for (auto& l : out)
    if (!l.starts_with("#")) l = l.substr(0, l.rfind(','));
  return out;
}
std::vector<std::string> without_notice(std::vector<std::string> rows) {
  std::erase_if(rows, [](const std::string& l) { return l.find("APPROXIMATE RESULT") != std::string::npos; });
  return rows;
}

// Solves both inputs with the same options on every available backend and
// compares the outputs. `notice_differs`: the Fluent case runs under a named
// approximation that the Kutta case, stating the setting, does not need; the
// outputs then differ only in the provenance notice.
void check_same_outputs(const std::filesystem::path& fluent, const std::string& fluent_options, const std::filesystem::path& kutta,
                        const std::string& kutta_options, bool notice_differs = false) {
  for (const auto backend : Device::available()) {
    const std::string name = backend_name(backend);
    INFO("backend: " << name);
    MeshTemp t;
    const auto common = " --iterations 3 --backend " + name + " ";
    REQUIRE_MESSAGE(kutta_cli("solve " + quote(fluent) + common + fluent_options + " --out " + quote(t.path("fluent")), t.path("fluent.log")) == 0,
                    mesh_text(t.path("fluent.log")));
    REQUIRE_MESSAGE(kutta_cli("solve " + quote(kutta) + common + kutta_options + " --out " + quote(t.path("kutta")), t.path("kutta.log")) == 0,
                    mesh_text(t.path("kutta.log")));
    for (const auto* file : {"fluxes.csv", "forces.csv", "result.vtu"}) {
      INFO(file);
      if (notice_differs)
        CHECK(without_notice(lines(t.path("fluent") / file)) == without_notice(lines(t.path("kutta") / file)));
      else
        CHECK(mesh_text(t.path("fluent") / file) == mesh_text(t.path("kutta") / file));
    }
    if (notice_differs) {
      auto a = std::get<Json::Object>(parse_json(mesh_text(t.path("fluent/loads.json"))).value);
      auto b = std::get<Json::Object>(parse_json(mesh_text(t.path("kutta/loads.json"))).value);
      CHECK(a.at("approximate") == Json{true});
      CHECK(b.at("approximate") == Json{false});
      for (auto* m : {&a, &b}) m->erase("approximate"), m->erase("provenance");
      CHECK(a == b);
      CHECK(without_notice(residual_rows(t.path("fluent/residuals.csv"))) == without_notice(residual_rows(t.path("kutta/residuals.csv"))));
    } else {
      CHECK(mesh_text(t.path("fluent/loads.json")) == mesh_text(t.path("kutta/loads.json")));
      CHECK(residual_rows(t.path("fluent/residuals.csv")) == residual_rows(t.path("kutta/residuals.csv")));
    }
    const auto summary = std::get<Json::Object>(parse_json(mesh_text(t.path("kutta/summary.json"))).value);
    CHECK(summary.at("case_format") == Json{std::string("kutta")});
  }
}

struct Pair {
  const char* fixture;
  const char* fluent_options;
  const char* kutta_options;
  bool notice_differs;
};
// The approximations are the Fluent case's. The Kutta case keeps the physics
// that needs them (WFGC on a polyhedral mesh, k-epsilon compressibility and
// Fluent's k-epsilon solver options), except the cylinder's: a Fluent 18 case
// leaves prevent-reverse-flow unwritten, the Kutta case states it.
constexpr Pair corpus_pairs[] = {{"mixing_elbow_v6", "--approximate wfgc", "--approximate wfgc", false},
                                 {"elbow_param", "--approximate compressibility,ke-options", "--approximate compressibility,ke-options", false},
                                 {"nozzle", "--approximate wfgc", "--approximate wfgc", false},
                                 {"cylinder_tutorial", "--approximate reverse-flow", "", true}};
setup::KuttaCase committed_case(const char* fixture) {
  return setup::read_kutta_case_file(std::filesystem::path(KUTTA_SOURCE_DIR "/tests/data/kutta_cases") / (std::string(fixture) + ".kutta.json"));
}
}  // namespace

TEST_CASE("The committed Kutta cases are what kutta convert writes for their Fluent cases") {
  for (const auto& p : corpus_pairs) {
    INFO(p.fixture);
    const auto source = committed_case(p.fixture);
    const auto fluent = fluent_examples() / source.mesh;
    if (!std::filesystem::exists(fluent)) {
      MESSAGE("SKIP: optional Fluent case absent: " << fluent.string());
      continue;
    }
    MeshTemp t;
    REQUIRE(kutta_cli("convert " + quote(fluent) + " --out " + quote(t.path("converted.kutta.json")), t.path("convert.log")) == 0);
    auto converted = setup::read_kutta_case_file(t.path("converted.kutta.json"));
    CHECK(std::filesystem::weakly_canonical(t.directory / converted.mesh) == std::filesystem::weakly_canonical(fluent));
    converted.mesh = source.mesh;
    CHECK(converted == source);
  }
}

TEST_CASE("Kutta cases converted from Fluent cases solve to the same bytes on every backend" * doctest::test_suite("slow")) {
  for (const auto& p : corpus_pairs) {
    INFO(p.fixture);
    auto source = committed_case(p.fixture);
    const auto fluent = fluent_examples() / source.mesh;
    if (!std::filesystem::exists(fluent)) {
      MESSAGE("SKIP: optional Fluent case absent: " << fluent.string());
      continue;
    }
    MeshTemp t;
    source.mesh = std::filesystem::absolute(fluent);
    mesh_text(t.path("case.kutta.json"), setup::kutta_case_text(source));
    check_same_outputs(fluent, p.fluent_options, t.path("case.kutta.json"), p.kutta_options, p.notice_differs);
  }
}

TEST_CASE("The NACA ground-effect wing from a bare mesh and a Kutta case solves like its donor Fluent case" * doctest::test_suite("slow")) {
  const auto donor = generated_fixture("wing/naca4412_ground.cas");
  const auto mesh = generated_fixture("wing/naca4412_ground.msh");
  const auto kutta = generated_fixture("wing/naca4412_ground.kutta.json");
  if (donor.empty() || mesh.empty() || kutta.empty()) return;
  check_same_outputs(donor, "", kutta, "");
}

TEST_CASE("kutta new-case writes a starter case for a bare mesh that solves") {
  MeshTemp t;
  std::filesystem::copy_file(KUTTA_SOURCE_DIR "/tests/data/meshes/fsae_box.msh", t.path("box.msh"));
  REQUIRE(kutta_cli("new-case " + quote(t.path("box.msh")) + " --speed 20 --role top=wall --out " + quote(t.path("box.kutta.json")) + " --json",
                    t.path("new.json")) == 0);
  const auto report = std::get<Json::Object>(parse_json(mesh_text(t.path("new.json"))).value);
  CHECK(report.at("speed") == Json{20.0});
  for (const auto& zone : std::get<Json::Array>(report.at("zones").value)) {
    const auto& z = std::get<Json::Object>(zone.value);
    if (z.at("name") == Json{std::string("top")}) {
      CHECK(z.at("role") == Json{std::string("wall")});
      CHECK(z.at("suggested_role") == Json{std::string("symmetry")});
    }
  }
  const auto c = setup::read_kutta_case_file(t.path("box.kutta.json"));
  CHECK(c.mesh == "box.msh");
  CHECK(std::holds_alternative<setup::WallInput>(c.boundaries.at("top")));
  // An existing case is not replaced without --force.
  CHECK(kutta_cli("new-case " + quote(t.path("box.msh")) + " --out " + quote(t.path("box.kutta.json")), t.path("again.log")) != 0);

  REQUIRE_MESSAGE(kutta_cli("solve " + quote(t.path("box.kutta.json")) + " --iterations 3 --backend cpu --out " + quote(t.path("run")),
                            t.path("solve.log")) == 0,
                  mesh_text(t.path("solve.log")));
  const auto summary = std::get<Json::Object>(parse_json(mesh_text(t.path("run/summary.json"))).value);
  CHECK(summary.at("case_format") == Json{std::string("kutta")});
  CHECK(summary.at("mesh_sha256") == Json{sha256_file(t.path("box.msh"))});
  const auto forces = mesh_text(t.path("run/forces.csv"));
  CHECK(forces.find("report,\"drag\"") != std::string::npos);
  CHECK(forces.find("report,\"downforce\"") != std::string::npos);
  // A Kutta case reads no Fluent data file.
  CHECK(kutta_cli("solve " + quote(t.path("box.kutta.json")) + " --data " + quote(t.path("box.msh")) + " --iterations 1 --out " +
                      quote(t.path("bad")),
                  t.path("bad.log")) != 0);
}

namespace {
// Independent native Fluent encodings: start with ordinary default settings,
// edit their Scheme records as Fluent's writer does, and save a real .cas.
// The test then exercises read_case -> convert CLI -> read JSON -> solve CLI.
void replace_entries(fluent::SExpr& alist, const std::string& text) {
  for (auto entry : fluent::parse_sexprs(text)) {
    const auto key = entry.elements.front().atom_text();
    const auto it = std::find_if(alist.elements.begin(), alist.elements.end(), [&](const auto& e) {
      return !e.elements.empty() && e.elements.front().atom_text() == key;
    });
    if (it == alist.elements.end()) alist.elements.push_back(std::move(entry));
    else *it = std::move(entry);
  }
}
void edit_thread(fluent::FluentCase& c, int id, const std::string& type, const std::string& entries) {
  auto& records = c.raw_sections.at("Thread Variables");
  const auto it = std::find_if(records.begin(), records.end(), [&](const auto& e) { return e.elements[1].elements[0].integer == id; });
  REQUIRE(it != records.end());
  it->elements[1].elements[1] = fluent::SExpr::symbol(type);
  replace_entries(it->elements[2], entries);
}
void write_native_case(const fluent::FluentCase& c, const std::filesystem::path& out) {
  auto text = mesh_text(KUTTA_SOURCE_DIR "/tests/data/meshes/fsae_box.msh");
  // Mesh section 13 stores hexadecimal zone IDs and BC type codes; section
  // 39 uses decimal IDs. Keep topology and Thread Variables consistent.
  for (const auto& e : c.raw_sections.at("Thread Variables")) {
    const auto id = e.elements[1].elements[0].integer;
    const auto type = *e.elements[1].elements[1].atom_text();
    if (type == "fluid") continue;
    const int code = type == "interior" ? 2 : type == "wall" ? 3 : type == "pressure-inlet" ? 4 : type == "pressure-outlet" ? 5 :
                     type == "symmetry" ? 7 : type == "velocity-inlet" ? 10 : type == "porous-jump" ? 14 : 20;
    std::ostringstream zone, hex_code;
    zone << std::hex << id;
    hex_code << std::hex << code;
    const std::regex header("(\\(13 \\(" + zone.str() + " [0-9a-f]+ [0-9a-f]+ )[0-9a-f]+( [0-9a-f]+\\))");
    std::smatch match;
    REQUIRE(std::regex_search(text, match, header));
    text.replace(std::size_t(match.position()), std::size_t(match.length()), match.str(1) + hex_code.str() + match.str(2));
  }
  const auto& threads = c.raw_sections.at("Thread Variables");
  if (std::any_of(threads.begin(), threads.end(), [](const auto& e) { return e.elements[1].elements[0].integer == 20; })) {
    const std::regex cells("\\(12 \\(2 1 46 1 0\\)\\([\\s\\S]*?\\)\\)");
    REQUIRE(std::regex_search(text, cells));
    text = std::regex_replace(text, cells, "(12 (2 1 18 1 4))\n(12 (14 19 46 1 4))");
  }
  text = std::regex_replace(text, std::regex("\\(45 [^\\n]+\\n"), "");
  text += "\n(0 \"fluent22.2\")\n";
  for (const auto& [section, records] : c.raw_sections)
    for (const auto& e : records) text += fluent::to_string(e) + '\n';
  mesh_text(out, text);
}
}

TEST_CASE("Native FSAE Fluent physics convert and solve to byte-identical Kutta outputs" * doctest::test_suite("slow")) {
  const auto mesh = read_mesh(KUTTA_SOURCE_DIR "/tests/data/meshes/fsae_box.msh");
  const auto info = setup::mesh_info(mesh);
  for (const auto* feature : {"rotating-wheel", "moving-components", "moving-direction", "pressure-normal", "pressure-direction",
                             "mass-rate", "mass-flux", "porous", "mrf", "relative-porous-mrf", "porous-jump", "transition-constants", "ke-options", "named-materials", "nested-mrf", "nested-mrf-absolute", "absolute-porous-mrf", "pressure-reference"}) {
    INFO(std::string(feature));
    MeshTemp t;
    auto input = setup::fsae_starter(info, "box.msh", {});
    if (std::string_view(feature).find("mrf") != std::string_view::npos)
      for (const auto* name : {"inlet", "outlet", "top", "side_left", "side_right"}) input.boundaries[name] = setup::WallInput{};
    if (std::string_view(feature) == "transition-constants") input.models.intermittency_transition = true;
    if (std::string_view(feature) == "ke-options") {
      input.models.viscous = setup::ViscousModel::k_epsilon_standard;
      input.methods.specific_dissipation_rate.reset();
      input.methods.turbulent_dissipation_rate = setup::Upwind::second_order;
    }
    auto native = setup::lower_kutta_case(input, info);
    native.raw_sections.at("Thread Variables").push_back(fluent::parse_sexpr("(39 (3 interior default-interior 1) ())"));
    const std::string name(feature);
    if (name == "rotating-wheel") edit_thread(native, 10, "wall", "(motion-bc . 1) (rotating? . #t) (relative? . #f) (omega . 12) (x-origin . .2) (y-origin . .1) (z-origin . .2) (ai . 0) (aj . -2) (ak . 0)");
    if (name == "moving-components") edit_thread(native, 6, "wall", "(motion-bc . 1) (rotating? . #f) (components? . #t) (relative? . #f) (u . 15) (v . 0) (w . 0)");
    if (name == "moving-direction") edit_thread(native, 6, "wall", "(motion-bc . 1) (rotating? . #f) (components? . #f) (relative? . #t) (vmag . 15) (ni . 3) (nj . 0) (nk . 0)");
    if (name.starts_with("pressure")) edit_thread(native, 4, "pressure-inlet", "(p0 . 150) (p . 2) (direction-spec . " + std::string(name == "pressure-normal" ? "1" : "0") + ") (ni . 2) (nj . 0) (nk . 0)");
    if (name.starts_with("mass")) edit_thread(native, 4, "mass-flow-inlet", name == "mass-rate" ? "(flow-spec . 0) (mass-flow . 3) (direction-spec . 1)" : "(flow-spec . 1) (mass-flux . 5) (direction-spec . 0) (ni . 2) (nj . 0) (nk . 0)");
    if (name == "porous" || name == "relative-porous-mrf" || name == "absolute-porous-mrf") edit_thread(native, 2, "fluid", R"((porous? . #t) (porosity . .35) (rel-vel-resistance? . #t)
      (porous-r/1 . 30) (porous-r/2 . 200) (porous-r/3 . 400) (porous-c/1 . 4) (porous-c/2 . 20) (porous-c/3 . 40)
      (direction-1/x . 1) (direction-1/y . 1) (direction-1/z . 0) (direction-2/x . -1) (direction-2/y . 1) (direction-2/z . 0)
      (alt-inertial-form? . #f) (c0 . 0) (c1 . 0))");
    if (name == "absolute-porous-mrf") edit_thread(native, 2, "fluid", "(rel-vel-resistance? . #f)");
    if (name.find("mrf") != std::string::npos) edit_thread(native, 2, "fluid", "(mrf-motion? . #t) (mrf-omega . .4) (mrf-grid-x-vel . .25) (mrf-x-origin . .2)");
    if (name == "porous-jump") edit_thread(native, 3, "porous-jump", "(alpha . .1) (dm . .02) (c2 . 1)");
    auto& rp = native.raw_sections.at("Rampant Variables")[0].elements[1];
    if (name == "pressure-reference") replace_entries(rp, "(pressure-reference/x . .3) (pressure-reference/y . .2) (pressure-reference/z . .1)");
    if (name == "transition-constants") replace_entries(rp, "(trans-one-eqn-ctu1 . 90) (trans-one-eqn-cpg2 . -.7) (intermit-sig . 1.2) (sst-limit-cross-diff-to-n-times-w-destruct . 123)");
    if (name == "ke-options") replace_entries(rp, "(ke-source-mod . #f) (ke-time-marching? . #f)");
    if (name == "named-materials" || name.starts_with("nested-mrf")) {
      auto& threads = native.raw_sections.at("Thread Variables");
      const auto base = std::find_if(threads.begin(), threads.end(), [](const auto& e) { return e.elements[1].elements[0].integer == 2; });
      REQUIRE(base != threads.end());
      auto child = *base;
      child.elements[1].elements[0] = fluent::parse_sexpr("20");
      child.elements[1].elements[2] = fluent::SExpr::symbol("fluid-copy");
      threads.push_back(std::move(child));
      if (name == "named-materials") {
        replace_entries(rp, "(materials ((air fluid (density . 1.225) (viscosity . 1.7894e-5)) (air-copy fluid (density . 1.225) (viscosity . 1.7894e-5))))");
        edit_thread(native, 20, "fluid", "(material . air-copy)");
      } else {
        edit_thread(native, 20, "fluid", "(mrf-relative-to-thread . 2) (mrf-omega . .1) (mrf-x-origin . .3) (mrf-grid-x-vel . .1)");
        if (name == "nested-mrf-absolute") replace_entries(rp, "(reference-frame/relative/origin-absolute? . #t)");
      }
    }
    write_native_case(native, t.path("native.cas"));
    REQUIRE_MESSAGE(kutta_cli("convert " + quote(t.path("native.cas")) + " --out " + quote(t.path("converted.kutta.json")), t.path("convert.log")) == 0,
                    mesh_text(t.path("convert.log")));
    check_same_outputs(t.path("native.cas"), "", t.path("converted.kutta.json"), "");
    std::cout << "KUTTA_NATIVE_EQUIVALENT " << feature << '\n';
  }
}

TEST_CASE("The FSAE wheel starter CLI saves geometry and solves on every backend") {
  MeshTemp t;
  std::filesystem::copy_file(KUTTA_SOURCE_DIR "/tests/data/meshes/fsae_wheel.msh", t.path("wheels.msh"));
  const auto args = "new-case " + quote(t.path("wheels.msh")) + " --speed 20 --wheel wheel_front=0.2,0.4,0.25,0,0,0,1 --out " + quote(t.path("wheels.kutta.json"));
  REQUIRE_MESSAGE(kutta_cli(args, t.path("new.log")) == 0, mesh_text(t.path("new.log")));
  const auto c = setup::read_kutta_case_file(t.path("wheels.kutta.json"));
  CHECK(std::get<setup::WallInput>(c.boundaries.at("wheel_front")).angular_speed == 100);
  for (const auto backend : Device::available()) {
    const auto name = std::string(backend_name(backend));
    REQUIRE_MESSAGE(kutta_cli("solve " + quote(t.path("wheels.kutta.json")) + " --iterations 3 --backend " + name + " --out " + quote(t.path(name)), t.path(name + ".log")) == 0,
                    mesh_text(t.path(name + ".log")));
    CHECK(std::filesystem::exists(t.path(name) / "loads.json"));
  }
  CHECK(kutta_cli("new-case " + quote(t.path("wheels.msh")) + " --wheel wheel_front=0,0,0,0,0,1,0", t.path("bad.log")) != 0);
  CHECK(kutta_cli("new-case " + quote(t.path("wheels.msh")) + " --wheel inlet=1,0,0,0,0,1,0", t.path("bad.log")) != 0);
  CHECK(kutta_cli("new-case " + quote(t.path("wheels.msh")) + " --wheel wheel_front=1,0,0,0,0,1,0,", t.path("bad.log")) != 0);
}

TEST_CASE("Corpus FSAE physics preserve strict solver holds or convert with identical outputs" * doctest::test_suite("slow")) {
  const auto root = fluent_examples();
  if (!std::filesystem::exists(root)) { MESSAGE("SKIP: optional Fluent corpus absent"); return; }
  std::vector<std::filesystem::path> cases;
  for (const auto& entry : std::filesystem::recursive_directory_iterator(root)) {
    const auto name = entry.path().filename().string();
    if (entry.is_regular_file() && (name.ends_with(".cas") || name.ends_with(".cas.h5"))) cases.push_back(entry.path());
  }
  std::sort(cases.begin(), cases.end());
  const auto flags = setup::parse_approximations("wfgc,compressibility,ke-options,reverse-flow");
  for (const auto& file : cases) {
    fluent::FluentCase c;
    try { c = fluent::read_case(file); }
    catch (const std::exception&) { continue; }  // malformed corpus files have separate reader gates
    const bool boundary = std::any_of(c.boundaries.begin(), c.boundaries.end(), [](const auto& z) {
      return z.type == fluent::ZoneType::pressure_inlet || z.type == fluent::ZoneType::mass_flow_inlet ||
             z.type == fluent::ZoneType::porous_jump || z.wall_motion == fluent::WallMotion::rotating || z.wall_motion == fluent::WallMotion::translating;
    });
    const bool zone = std::any_of(c.cell_zones.begin(), c.cell_zones.end(), [](const auto& z) { return z.frame_motion == true || z.porous == true; });
    std::set<std::string> fluid_materials;
    for (const auto& z : c.cell_zones) if (z.type == fluent::ZoneType::fluid && z.material) fluid_materials.insert(*z.material);
    const auto stored = [&](const std::string& key) { return c.variables.value(key); };
    const auto false_flag = [&](const std::string& key) {
      const auto* v = stored(key);
      return v && v->kind == fluent::SExpr::Kind::boolean && !v->boolean;
    };
    const bool model = c.viscous.transition == fluent::TransitionModel::intermittency ||
      (c.viscous.model == fluent::ViscousModel::ke_standard && (false_flag("ke-source-mod") || false_flag("ke-time-marching?"))) ||
      (c.viscous.model == fluent::ViscousModel::kw_sst && stored("sst-limit-cross-diff-to-n-times-w-destruct") &&
       stored("sst-limit-cross-diff-to-n-times-w-destruct")->number() != 1e6);
    bool reference = false;
    for (const auto* axis : {"x", "y", "z"}) {
      const auto* v = stored(std::string("pressure-reference/") + axis);
      reference = reference || (v && v->number() && *v->number() != 0);
    }
    if (!boundary && !zone && fluid_materials.size() < 2 && !model && !reference) continue;
    const auto relative = file.lexically_relative(root).generic_string();
    INFO(relative);
    std::string options;
    try {
      const auto problem = setup::fluent_setup(c, read_mesh(file), flags);
      for (const auto flag : problem.approximations) options += (options.empty() ? "--approximate " : ",") + std::string(setup::approximation_info(flag).name);
    } catch (const std::exception& e) {
      // Do not turn an original .cas rejection into a claimed equivalence.
      // Keep the exact original blocker visible in slow-tier evidence.
      std::cout << "KUTTA_CORPUS_HOLD " << relative << '\n' << e.what() << '\n';
      continue;
    }
    MeshTemp t;
    REQUIRE_MESSAGE(kutta_cli("convert " + quote(file) + " --out " + quote(t.path("converted.kutta.json")), t.path("convert.log")) == 0,
                    mesh_text(t.path("convert.log")));
    check_same_outputs(file, options, t.path("converted.kutta.json"), options);
    std::cout << "KUTTA_CORPUS_EQUIVALENT " << relative << '\n';
  }
}

TEST_CASE("Native Fluent FSAE boundary topology codes must match thread types") {
  const auto info = setup::mesh_info(read_mesh(KUTTA_SOURCE_DIR "/tests/data/meshes/fsae_box.msh"));
  for (const auto* type : {"pressure-inlet", "mass-flow-inlet", "porous-jump"}) {
    INFO(std::string(type));
    MeshTemp t;
    auto native = setup::lower_kutta_case(setup::fsae_starter(info, "box.msh", {}), info);
    native.raw_sections.at("Thread Variables").push_back(fluent::parse_sexpr("(39 (3 interior default-interior 1) ())"));
    const bool jump = std::string_view(type) == "porous-jump";
    edit_thread(native, jump ? 3 : 4, type, jump ? "(alpha . .1) (dm . .02) (c2 . 1)" :
                std::string_view(type) == "pressure-inlet" ? "(p0 . 1) (direction-spec . 1)" : "(flow-spec . 0) (mass-flow . .1) (direction-spec . 1)");
    write_native_case(native, t.path("native.cas"));
    const auto issue_path = std::string("mesh/zone/") + (jump ? "3" : "4") + "/type";
    const auto valid = fluent::read_case(t.path("native.cas"));
    CHECK(std::none_of(valid.issues.begin(), valid.issues.end(), [&](const auto& i) { return i.path == issue_path; }));
    auto text = mesh_text(t.path("native.cas"));
    const std::regex header(std::string("(\\(13 \\(") + (jump ? "3" : "4") + " [0-9a-f]+ [0-9a-f]+ )[0-9a-f]+( [0-9a-f]+\\))");
    std::smatch match;
    REQUIRE(std::regex_search(text, match, header));
    text.replace(std::size_t(match.position()), std::size_t(match.length()), match.str(1) + "3" + match.str(2));
    mesh_text(t.path("bad.cas"), text);
    const auto invalid = fluent::read_case(t.path("bad.cas"));
    CHECK(std::any_of(invalid.issues.begin(), invalid.issues.end(), [&](const auto& i) { return i.path == issue_path; }));
  }
}
