// Kutta cases (mesh + *.kutta.json) through the real CLI: a case converted from
// a Fluent case solves to the same bytes as the Fluent case, and a starter
// written from a bare mesh runs.
#include "../backend_helpers.hpp"
#include "../fluent_examples.hpp"
#include "../generated_fixtures.hpp"
#include "../io/test_support.hpp"
#include "kutta/core/json_value.hpp"
#include "kutta/setup/kutta_case.hpp"
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
    REQUIRE(kutta_cli("solve " + quote(fluent) + common + fluent_options + " --out " + quote(t.path("fluent")), t.path("fluent.log")) == 0);
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
