#include "kutta/case/fluent_case.hpp"
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <doctest/doctest.h>
#include <fstream>
#include <tuple>

using namespace kutta::fluent;
namespace {
std::filesystem::path examples() {
  auto env = std::getenv("KUTTA_FLUENT_EXAMPLES");
  return env ? env : KUTTA_SOURCE_DIR "/validation/fluent_examples";
}
bool available(const std::filesystem::path &path) {
  if (std::filesystem::exists(path))
    return true;
  MESSAGE("Fluent fixture absent, optional real-file test skipped: ", path.string());
  return false;
}
const BoundaryZone &boundary(const FluentCase &c, std::int64_t id) {
  auto it =
      std::find_if(c.boundaries.begin(), c.boundaries.end(), [id](auto &z) { return z.id == id; });
  if (it == c.boundaries.end())
    throw std::runtime_error("missing expected boundary");
  return *it;
}
FluentCase synthetic(std::string rp, std::string threads = {}) {
  std::map<std::string, std::vector<SExpr>> sections;
  sections["Rampant Variables"] = parse_sexprs("(37 (" + rp + "))");
  sections["Thread Variables"] = parse_sexprs(threads);
  return extract_case(std::move(sections));
}
std::string config(bool ke, bool kw, bool lam = false) {
  return "(case-config ((rp-3d? . #t) (rp-double? . #t) (rp-unsteady? . #f) (rp-seg? . #t) "
         "(rp-lam? . " +
         std::string(lam ? "#t" : "#f") + ") (rp-ke? . " + (ke ? "#t" : "#f") + ") (rp-kw? . " +
         (kw ? "#t" : "#f") + "))) ";
}
auto typed(const Material &z) {
  return std::tie(z.name, z.type, z.density, z.viscosity, z.nonconstant_density,
                  z.nonconstant_viscosity);
}
auto typed(const BoundaryZone &z) {
  return std::tie(
      z.id, z.domain, z.name, z.fluent_type, z.type, z.velocity_spec, z.velocity_magnitude,
      z.gauge_pressure, z.velocity, z.direction, z.turbulence_spec, z.turbulence, z.wall_motion,
      z.relative_motion, z.specified_shear, z.prevent_reverse_flow, z.frame_of_reference,
      z.backflow_pressure_spec, z.wall_components, z.shear_bc, z.rough_bc, z.angular_velocity, z.roughness_height,
      z.roughness_constant, z.axis_origin, z.axis_direction, z.shear);
}
auto typed(const CellZone &z) {
  return std::tie(z.id, z.domain, z.name, z.fluent_type, z.type, z.material, z.frame_motion,
                  z.porous, z.laminar, z.relative_resistance, z.relative_to_zone, z.motion_spec,
                  z.angular_velocity, z.porosity, z.frame_velocity, z.axis_origin, z.axis_direction,
                  z.viscous_resistance, z.inertial_resistance, z.porous_directions);
}
auto typed(const ReportDefinition &r) {
  return std::tie(r.name, r.fluent_type, r.kind, r.zone_ids, r.surface_ids, r.reported_zone_ids, r.scaled, r.zone_names, r.direction);
}
template <class T> void same(const std::vector<T> &a, const std::vector<T> &b) {
  REQUIRE(a.size() == b.size());
  for (std::size_t i = 0; i < a.size(); ++i) {
    INFO("record ", i, " expected ", a[i].name, " actual ", b[i].name);
    CHECK(typed(a[i]) == typed(b[i]));
  }
}
void same_settings(const FluentCase &a, const FluentCase &b) {
  CHECK(a.solver_version == b.solver_version);
  CHECK(a.dimension == b.dimension);
  CHECK(a.double_precision == b.double_precision);
  CHECK(a.unsteady == b.unsteady);
  CHECK(a.pressure_based == b.pressure_based);
  CHECK(a.viscous == b.viscous);
  CHECK(a.numerics == b.numerics);
  CHECK(a.operating_pressure == b.operating_pressure);
  CHECK(a.reference == b.reference);
  same(a.materials, b.materials);
  same(a.boundaries, b.boundaries);
  same(a.cell_zones, b.cell_zones);
  same(a.reports, b.reports);
}
} // namespace
TEST_CASE("Case extraction uses authoritative flags and never invents missing settings") {
  auto c =
      synthetic(config(false, true) +
                "(rp-kw? #f) (kw-sst-on? #t) (kw-bsl-on? #f) (kw-std-on? #f) (flow/scheme 24) "
                "(flow/scheme 20 integer #f 20) (pressure/scheme 999) (unknown-setting '(a . b))");
  CHECK(c.viscous.model == ViscousModel::kw_sst);
  CHECK(c.numerics.coupling == Coupling::coupled);
  CHECK_FALSE(c.numerics.pressure);
  CHECK_FALSE(c.operating_pressure);
  CHECK_FALSE(c.reference.area);
  REQUIRE(c.variables.value("unknown-setting"));
  CHECK(to_string(*c.variables.value("unknown-setting")) == "(quote (a . b))");
  auto raw = c.variables.unparsed_names();
  CHECK(std::find(raw.begin(), raw.end(), "unknown-setting") != raw.end());
  CHECK(std::any_of(c.issues.begin(), c.issues.end(), [](auto &i) {
    return i.path == "rp/pressure/scheme" && i.reason.starts_with("unparsed");
  }));
  CHECK(synthetic(config(false, false, true)).viscous.model == ViscousModel::laminar);
  CHECK(synthetic(config(true, false) + "(rng-ke-on? #t) (realizable-epsilon? #f)").viscous.model ==
        ViscousModel::ke_rng);
  CHECK(synthetic(config(true, false) + "(rng-ke-on? #f) (realizable-epsilon? #t)").viscous.model ==
        ViscousModel::ke_realizable);
  CHECK(synthetic(config(true, false) +
                  "(rng-ke-on? #f) (realizable-epsilon? #f) (ke-realizability-on? #t)")
            .viscous.model == ViscousModel::ke_standard);
  CHECK(synthetic(config(false, true) + "(kw-sst-on? #f) (kw-bsl-on? #t) (kw-std-on? #f)")
            .viscous.model == ViscousModel::kw_bsl);
  CHECK(synthetic(config(false, true) + "(kw-sst-on? #f) (kw-bsl-on? #f) (kw-std-on? #t)")
            .viscous.model == ViscousModel::kw_standard);
  CHECK_FALSE(synthetic(config(true, true)).viscous.model);
}
TEST_CASE("Case expression selectors do not expose inactive numeric settings") {
  auto c = synthetic(config(false, false, true) + R"scm(
    (operating-pressure 101325) (operating-pressure-sel "expression")
    (operating-pressure-expr "p(t)")
    (pseudo-auto-time-step-scale-factor 1)
    (pseudo-auto-time-step-scale-factor-sel "expression")
    (pseudo-auto-time-step-scale-factor-expr "scale(t)")
  )scm");
  CHECK_FALSE(c.operating_pressure);
  CHECK_FALSE(c.numerics.pseudo.scale_factor);
  REQUIRE(c.variables.value("operating-pressure-expr"));
  CHECK(c.variables.value("operating-pressure-expr")->text == "p(t)");
}
TEST_CASE("Case BC material frame porous and force bindings stay explicit") {
  auto c = synthetic(config(false, true) + R"scm(
    (kw-sst-on? #t) (kw-bsl-on? #f) (kw-std-on? #f)
    (kw-low-re-mod? #t) (kw-sst-strain? #f) (kw-not-trans-sst-production-limiter? #t)
    (sst-pk-factor 10) (turbulence/wall-omega-treatment tabulated)
    (materials ((water fluid (density (constant . 998)) (viscosity (sutherland 1 2 3) (constant . 1)))))
    (monitor/report-definitions (((name . "drag") (report-definition drag (zone-list 8 9) (force-vector 1 0 0)))))
  )scm",
                     R"scm(
    (39 (4 velocity-inlet inlet 1) ((velocity-spec . 2) (vmag (expression . "speed(t)") (constant . 99)) (ke-spec . 2) (turb-intensity . 0.05) (turb-viscosity-ratio . 10)))
    (39 (8 wall wheel 1) ((motion-bc . 1) (rotating? . #t) (omega . 100) (relative? . #f) (ai . 0) (aj . 0) (ak . 1) (specified-shear? . #t) (shear-x . 2) (roughness-height . 0.001) (roughness-const . 0.5)))
    (39 (9 fluid rotor 1) ((mrf-motion? . #t) (mrf-omega (constant . 100) (profile "" "")) (porous? . #t) (porous-r/1 (constant . 200)) (porous-c/1 . 4) (direction-1/x . 1) (reference-frame-axis-direction-components ((constant . 0)) ((constant . 0)) ((constant . 1)))))
  )scm");
  CHECK(c.viscous.wall_omega == WallOmega::tabulated);
  CHECK(c.viscous.low_re == true);
  CHECK(c.viscous.production_factor == 10);
  auto &inlet = boundary(c, 4);
  REQUIRE(inlet.velocity_magnitude);
  CHECK(inlet.velocity_magnitude->kind == ScalarValue::Kind::expression);
  CHECK_FALSE(inlet.velocity_magnitude->constant);
  CHECK(inlet.velocity_magnitude->binding == "speed(t)");
  auto &wall = boundary(c, 8);
  CHECK(wall.wall_motion == WallMotion::rotating);
  CHECK(wall.angular_velocity->constant == 100);
  CHECK(wall.specified_shear == true);
  CHECK(wall.roughness_height->constant == 0.001);
  CHECK(wall.shear[0]->constant == 2);
  REQUIRE(c.cell_zones.size() == 1);
  auto &zone = c.cell_zones[0];
  CHECK(zone.frame_motion == true);
  CHECK(zone.angular_velocity->constant == 100);
  CHECK(zone.viscous_resistance[0]->constant == 200);
  CHECK(zone.inertial_resistance[0]->constant == 4);
  CHECK(zone.axis_direction[2]->constant == 1);
  REQUIRE(c.reports.size() == 1);
  CHECK(c.reports[0].kind == ReportDefinition::Kind::drag);
  CHECK_FALSE(c.reports[0].zone_ids);
  REQUIRE(c.reports[0].reported_zone_ids);
  CHECK(*c.reports[0].reported_zone_ids == std::vector<std::int64_t>{8, 9});
  REQUIRE(c.reports[0].direction);
  CHECK((*c.reports[0].direction)[0] == 1);
  REQUIRE(c.materials.size() == 1);
  CHECK(c.materials[0].nonconstant_viscosity == true);
  CHECK_FALSE(c.materials[0].viscosity->constant);
}
TEST_CASE("Case malformed identities and incomplete Scheme fail at the boundary") {
  CHECK_THROWS(synthetic(config(false, true), "(39 (bad wall x 1) ())"));
  CHECK_THROWS(synthetic(config(false, true), "(39 (8 wall a 1) ()) (39 (8 wall b 1) ())"));
  CHECK_THROWS(extract_case({}));
  CHECK_THROWS(synthetic("(x . )"));
}
TEST_CASE("Case legacy metadata is hexadecimal and binary bodies stay opaque") {
  auto path =
      std::filesystem::temp_directory_path() /
      ("kutta-e-legacy-" +
       std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()) + ".cas");
  struct Cleanup {
    std::filesystem::path path;
    ~Cleanup() {
      std::error_code e;
      std::filesystem::remove(path, e);
    }
  } cleanup{path};
  auto write = [&](std::string body) {
    std::ofstream f(path, std::ios::binary);
    f << "(0 \"fluent24.1.0 build\") (2 3) (37 (" << config(false, false, true) << ")) " << body;
  };
  write("(12 (0 1 10 0)) (12 (a 1 10 1 4)) (13 (0 1 20 0)) (3013 (b 1 20 3 4) opaque body End of "
        "Binary Section   3013) (39 (10 fluid f 1) ()) (39 (11 wall w 1) ())");
  auto c = read_case(path);
  CHECK(c.mesh.cells == 16);
  CHECK(c.mesh.faces == 32);
  CHECK(c.mesh.cell_zones[0].zone_id == 10);
  CHECK(c.mesh.face_zones[0].name == "w");
  write("(3013 (b 1 20 3 4) opaque End of Binary Section 3010)");
  CHECK_THROWS(read_case(path));
  write("(13 (0 1 20 0)) (13 (b 1 10 3 4))");
  CHECK_THROWS(read_case(path));
}
TEST_CASE("Real legacy cases parse without HDF5") {
  auto root = examples() / "static_mixer";
  for (auto n :
       {"Static_Mixer_Parameters_legacy_text.cas", "Static_Mixer_Parameters_legacy_binary.cas"}) {
    auto path = root / n;
    if (!available(path))
      return;
    auto c = read_case(path);
    CHECK(c.viscous.model == ViscousModel::kw_sst);
    CHECK(c.numerics.coupling == Coupling::coupled);
    CHECK(boundary(c, 17).velocity_magnitude->resolved == 1);
    CHECK(c.reports.size() == 2);
  }
}
#if KUTTA_HAS_HDF5
TEST_CASE("Real mixing elbow typed settings match verified CFF values") {
  auto path = examples() / "mixing_elbow/mixing_elbow.cas.h5";
  if (!available(path))
    return;
  auto c = read_case(path);
  REQUIRE(c.energy == true);
  CHECK(c.variables.value("rf-energy?")->boolean == false);
  CHECK(std::any_of(c.issues.begin(), c.issues.end(), [](auto &i) {
    return i.path == "case-config/rf-energy?" && i.classification == IssueClass::informational &&
           i.reason.starts_with("energy ignored (decoupled)");
  }));
  CHECK(c.solver_version == "22.2");
  CHECK(c.dimension == 3);
  CHECK(c.viscous.model == ViscousModel::kw_sst);
  CHECK(c.mesh.cells == 17822);
  CHECK(c.mesh.faces == 91581);
  CHECK(c.mesh.nodes == 66417);
  CHECK(c.numerics.coupling == Coupling::coupled);
  CHECK(c.numerics.pressure == Scheme::second_pressure);
  CHECK(c.numerics.k == Scheme::second_upwind);
  CHECK(c.numerics.omega == Scheme::second_upwind);
  CHECK(boundary(c, 30).name == "hot-inlet");
  CHECK(boundary(c, 30).type == ZoneType::velocity_inlet);
  CHECK(boundary(c, 30).velocity_spec == VelocitySpec::normal);
  CHECK(boundary(c, 30).velocity_magnitude->constant == 1.2);
  CHECK(boundary(c, 31).name == "cold-inlet");
  CHECK(boundary(c, 31).velocity_magnitude->constant == 1);
  CHECK(boundary(c, 30).turbulence.at("turb-hydraulic-diam").kind == ScalarValue::Kind::expression);
  CHECK(boundary(c, 30).turbulence.at("turb-hydraulic-diam").binding == "1 [in]");
  CHECK_FALSE(boundary(c, 30).turbulence.at("turb-hydraulic-diam").constant);
  CHECK(boundary(c, 30).turbulence.at("turb-hydraulic-diam").resolved == doctest::Approx(0.0254));
  CHECK(boundary(c, 32).gauge_pressure->constant == 0);
  CHECK(c.operating_pressure == 101325);
  CHECK(c.reference.density == 1.225);
  CHECK(c.reference.viscosity == 1.7894e-5);
  CHECK_FALSE(c.viscous.wall_omega);
  CHECK(c.cell_zones[0].relative_to_zone == -1);
  CHECK(std::none_of(c.issues.begin(), c.issues.end(), [](auto &i) {
    return i.path == "zone/87/mrf-relative-to-thread";
  }));
}
TEST_CASE("Real elbow1 authoritative k epsilon and velocity components") {
  auto path = examples() / "file_session/elbow1.cas.h5";
  if (!available(path))
    return;
  auto c = read_case(path);
  CHECK(c.solver_version == "24.1");
  CHECK(c.viscous.model == ViscousModel::ke_standard);
  REQUIRE(c.variables.value("rp-ke?"));
  CHECK_FALSE(c.variables.value("rp-ke?")->boolean);
  CHECK(c.config.value("rp-ke?")->boolean);
  CHECK(c.viscous.wall_omega == WallOmega::correlation);
  CHECK(c.mesh.cells == 13852);
  CHECK(c.mesh.faces == 44500);
  CHECK(c.numerics.coupling == Coupling::simple);
  CHECK(c.numerics.momentum == Scheme::first_upwind);
  CHECK(c.numerics.epsilon == Scheme::first_upwind);
  CHECK_FALSE(c.numerics.omega);
  CHECK(c.numerics.relaxation.contains("epsilon"));
  CHECK(c.viscous.unsupported_options.at("ke-std-wall?") == true);
  CHECK(c.viscous.unsupported_options.at("ke-enh-wall?") == false);
  CHECK(boundary(c, 6).velocity[1]->constant == 1.2);
  CHECK(*boundary(c, 7).velocity[0]->constant == doctest::Approx(0.40000001));
  CHECK(c.cell_zones[0].material == "water");
}
TEST_CASE("Static mixer typed settings agree across CFF ASCII and binary legacy") {
  auto root = examples() / "static_mixer";
  for (auto n : {"Static_Mixer_Parameters.cas.h5", "Static_Mixer_Parameters_legacy_text.cas",
                 "Static_Mixer_Parameters_legacy_binary.cas"})
    if (!available(root / n))
      return;
  auto c = read_case(root / "Static_Mixer_Parameters.cas.h5");
  auto a = read_case(root / "Static_Mixer_Parameters_legacy_text.cas");
  auto b = read_case(root / "Static_Mixer_Parameters_legacy_binary.cas");
  CHECK(c.viscous.model == ViscousModel::kw_sst);
  CHECK(c.mesh.cells == 22771);
  CHECK(c.mesh.faces == 119140);
  CHECK(c.reports.size() == 2);
  CHECK(boundary(c, 17).velocity_magnitude->kind == ScalarValue::Kind::parameter);
  CHECK(boundary(c, 17).velocity_magnitude->binding == "real-1");
  CHECK(boundary(c, 17).velocity_magnitude->resolved == 1);
  same_settings(c, a);
  same_settings(c, b);
  // Writer provenance differs; all typed physical/control values and ID intervals must match.
  CHECK(a.origin != c.origin);
  CHECK(a.mesh.cells == c.mesh.cells);
  CHECK(a.mesh.faces == c.mesh.faces);
  CHECK(a.mesh.nodes == c.mesh.nodes);
  CHECK(a.mesh.face_zones == c.mesh.face_zones);
  CHECK(b.mesh.face_zones == c.mesh.face_zones);
  CHECK(b.mesh.cell_zones == c.mesh.cell_zones);
  CHECK(a.raw_sections.at("Cortex Variables")[0].elements[0].integer == 38);
  REQUIRE(a.domains.contains(3));
  REQUIRE(a.domains.at(3).value("pb-qmom/n-pb-qmom"));
  CHECK(a.domains.at(3).value("pb-qmom/n-pb-qmom")->integer == 6);
}
#else
TEST_CASE("No HDF5 build reports unsupported files clearly") {
  CHECK_THROWS_WITH(read_case("unavailable.cas.h5"),
                    doctest::Contains("HDF5 support is not built"));
}
#endif

namespace {
bool issue_at(const FluentCase &c, std::string_view path) {
  return std::any_of(c.issues.begin(), c.issues.end(), [&](auto &i) { return i.path == path; });
}
}
TEST_CASE("Case support is an allowlist for authoritative physics flags") {
  for (auto key : {"rp-trans-sst?", "rp-sas?", "rp-des?", "rp-les?", "rp-axi?",
                   "rp-inviscid?", "sg-swirl?", "sg-mphase?", "sg-dynmesh?", "future-model?"}) {
    auto rp = config(false, true);
    rp.insert(rp.size() - 3, " (" + std::string(key) + " . #t)");
    auto c = synthetic(rp + "(kw-sst-on? #t) (kw-bsl-on? #f) (kw-std-on? #f)");
    INFO(key);
    CHECK(issue_at(c, "case-config/" + std::string(key)));
  }
  for (auto [key, val] : {std::pair{"rp-3d?", "#f"}, {"rp-visc?", "#f"},
                           {"rp-unsteady?", "#t"}, {"rp-seg?", "#f"}}) {
    auto c = synthetic("(case-config ((" + std::string(key) + " . " + val + "))) ");
    CHECK(issue_at(c, "case-config/" + std::string(key)));
  }
  auto c = synthetic(config(false, true) + "(kw-sst-on? #t) (kw-bsl-on? #f) (kw-std-on? #f) (ke-easm-on? #t) (gravity? #t) (gravity/x 0) (gravity/y -9.81) (gravity/z 0)");
  CHECK(issue_at(c, "rp/ke-easm-on?"));
  CHECK(issue_at(c, "rp/gravity?"));
}
TEST_CASE("Unsupported k epsilon retains its active equation schemes") {
  auto c = synthetic(config(true, false) + "(rng-ke-on? #f) (realizable-epsilon? #f) (epsilon/scheme 0) (epsilon/relax 0.8) (epsilon/pseudo-relax 0.75)");
  CHECK(issue_at(c, "viscous-model"));
  CHECK(c.variables.interpreted.contains("epsilon/scheme"));
  CHECK(c.variables.interpreted.contains("epsilon/relax"));
  CHECK(c.variables.interpreted.contains("epsilon/pseudo-relax"));
  CHECK_FALSE(issue_at(c, "rp/omega/scheme"));
}
TEST_CASE("Node reconstruction mapping remains explicitly unverified") {
  auto c = synthetic(config(false, false, true) + "(recon/cell-lsq? #f) (recon/node-lsq? #t)");
  CHECK(issue_at(c, "rp/recon/node-lsq?"));
}

TEST_CASE("Energy and gravity require proven constant flow properties") {
  for (bool constant : {false, true}) {
    auto rp = config(false, false, true);
    rp.insert(rp.size() - 3, " (rf-energy? . #t)");
    rp += "(gravity? #t) (gravity/x 0) (gravity/y -9.81) (gravity/z 0) (materials ((water fluid (density (";
    rp += constant ? "constant . 998" : "ideal-gas";
    rp += ")) (viscosity (constant . 0.001)))))";
    auto c = synthetic(rp, "(39 (2 fluid fluid 1) ((material . water)))");
    CHECK(c.gravity_vector[1] == -9.81);
    for (auto path : {"case-config/rf-energy?", "rp/gravity?"}) {
      auto i = std::find_if(c.issues.begin(), c.issues.end(), [&](auto &i) { return i.path == path; });
      REQUIRE(i != c.issues.end());
      CHECK(i->classification == (constant ? IssueClass::informational : IssueClass::unsupported));
    }
  }
}
TEST_CASE("Zone switches reject nonneutral cell inlet and outlet behavior") {
  for (auto key : {"sources?", "fixed?", "fanzone?", "mgrid-motion?", "solid-motion?", "deactivated?", "conical?", "future-switch?"}) {
    auto c = synthetic(config(false, false, true), "(39 (2 fluid f 1) ((" + std::string(key) + " . #t)))");
    CHECK(issue_at(c, "zone/2/" + std::string(key)));
  }
  for (auto [key, val] : {std::pair{"avg-press-spec?", "#t"}, {"targeted-mf-boundary?", "#t"},
      {"radial?", "#t"}, {"gen-nrbc-spec", "1"}, {"p-backflow-spec", "1"}}) {
    auto c = synthetic(config(false, false, true), "(39 (5 pressure-outlet out 1) ((" + std::string(key) + " . " + val + ")))");
    CHECK(issue_at(c, "zone/5/" + std::string(key)));
  }
  for (auto key : {"coordinate-system", "omega-swirl"}) {
    auto c = synthetic(config(false, false, true), "(39 (4 velocity-inlet in 1) ((" + std::string(key) + " . 1)))");
    CHECK(issue_at(c, "zone/4/" + std::string(key)));
  }
}
TEST_CASE("Translating walls require only the active velocity representation") {
  for (bool components : {false, true}) {
    auto settings = std::string("(motion-bc . 1) (rotating? . #f) (components? . ") + (components ? "#t" : "#f") + ") ";
    settings += components ? "(u . 15) (v . 0) (w . 0)" : "(vmag . 15) (ni . 1) (nj . 0) (nk . 0)";
    auto c = synthetic(config(false, false, true), "(39 (8 wall ground 1) (" + settings + "))");
    auto &z = boundary(c, 8);
    CHECK(z.wall_motion == WallMotion::translating);
    CHECK(z.wall_components == components);
    CHECK(z.velocity_spec == (components ? VelocitySpec::components : VelocitySpec::magnitude_direction));
    if (components) {
      REQUIRE(z.velocity[0]);
      CHECK(z.velocity[0]->constant == 15);
      CHECK_FALSE(z.velocity_magnitude);
    } else {
      REQUIRE(z.velocity_magnitude);
      CHECK(z.velocity_magnitude->constant == 15);
      CHECK_FALSE(z.velocity[0]);
    }
    for (auto key : {"u", "v", "w", "vmag", "ni", "nj", "nk"})
      CHECK_FALSE(issue_at(c, "zone/8/" + std::string(key)));
  }
}
TEST_CASE("Inlet singular vector override detects conflicting scalar components") {
  auto c = synthetic(config(false, false, true), R"scm(
    (39 (4 velocity-inlet in 1) ((velocity-spec . 1) (u . 2) (v . 0) (w . 0)
       (velocity-component ((constant . 3)) ((constant . 0)) ((constant . 0)))))
  )scm");
  CHECK(boundary(c, 4).velocity[0]->constant == 3);
  CHECK(issue_at(c, "zone/4/velocity-component"));
}
TEST_CASE("Parameter bindings never use their inactive saved constants") {
  auto c = synthetic(config(false, false, true) + R"scm(
    (named-expressions (((name . "speed") (parameterid . "real-1") (definition . "1 [m/s]"))
                        ((name . "dynamic") (parameterid . "real-2") (definition . "sin(t) [m/s]"))))
  )scm", R"scm(
    (39 (4 velocity-inlet in 1) ((velocity-spec . 2) (vmag (parameter . real-1) (constant . 2))))
    (39 (5 velocity-inlet other 1) ((velocity-spec . 2) (vmag (parameter . real-2) (constant . 99))))
  )scm");
  CHECK_FALSE(boundary(c, 4).velocity_magnitude->constant);
  CHECK(boundary(c, 4).velocity_magnitude->resolved == 1);
  CHECK(boundary(c, 4).velocity_magnitude->expression_name == "speed");
  CHECK_FALSE(boundary(c, 5).velocity_magnitude->constant);
  CHECK_FALSE(issue_at(c, "zone/4/vmag"));
  CHECK(issue_at(c, "zone/5/vmag"));
}
TEST_CASE("Only complete constant expressions with recognized units resolve") {
  for (auto expression : {"2 [m/s]", "2 [m]", "2 [Pa]", "2 [kg/m^3]", "2 [rad/s]", "2 [in]", "2 [mm]"}) {
    auto c = synthetic(config(false, false, true), "(39 (4 velocity-inlet in 1) ((velocity-spec . 2) (vmag (expression . \"" + std::string(expression) + "\"))))");
    CHECK_FALSE(boundary(c, 4).velocity_magnitude->constant);
    CHECK_FALSE(issue_at(c, "zone/4/vmag"));
    const auto expected = std::string_view(expression).ends_with("[in]") ? 0.0508
        : (std::string_view(expression).ends_with("[mm]") ? 0.002 : 2.0);
    CHECK(boundary(c, 4).velocity_magnitude->resolved == doctest::Approx(expected));
  }
  for (auto expression : {"sin(t)", "1 [furlong]", "1 [m] + 2 [m]", "1 [m] trailing", "1e999 [m]"}) {
    auto c = synthetic(config(false, false, true), "(39 (4 velocity-inlet in 1) ((velocity-spec . 2) (vmag (expression . \"" + std::string(expression) + "\"))))");
    CHECK(issue_at(c, "zone/4/vmag"));
  }
}
TEST_CASE("Force report surface IDs never masquerade as thread IDs") {
  auto c = synthetic(config(false, false, true) + R"scm(
    (monitor/report-definitions (
      ((name . "drag") (report-definition force-drag (surface-ids 1) (surface-names body) (force-vector 1 0 0) (scaled? . #t)))
      ((name . "force") (report-definition force-force (surface-ids 1) (force-vector 0 1 0) (scaled? . #f)))
      ((name . "bad") (report-definition force-lift (surface-names missing) (force-vector 0 1 0) (scaled? . #t)))))
  )scm", "(39 (8 wall body 1) ())");
  const auto &drag = c.reports[1];
  REQUIRE(drag.zone_ids);
  CHECK(*drag.zone_ids == std::vector<std::int64_t>{8});
  CHECK_FALSE(issue_at(c, "report/drag/zones"));
  CHECK_FALSE(c.reports[2].zone_ids);
  CHECK(issue_at(c, "report/force/zones"));
  CHECK(issue_at(c, "report/bad/zones"));
  CHECK(drag.raw.interpreted.contains("scaled?"));
  CHECK(drag.scaled == true);
  CHECK(c.reports[2].scaled == false);
  REQUIRE(drag.surface_ids);
  CHECK(*drag.surface_ids == std::vector<std::int64_t>{1});
}

TEST_CASE("Inactive turbulence flags remain visible in a laminar case") {
  auto c = synthetic(config(false, false, true) + "(kw-future? #t)");
  auto i = std::find_if(c.issues.begin(), c.issues.end(), [](auto &i) { return i.path == "rp/kw-future?"; });
  REQUIRE(i != c.issues.end());
  CHECK(i->classification == IssueClass::informational);
}
TEST_CASE("Active pseudo relaxation requires turbulence factors") {
  auto c = synthetic(config(false, true) + "(kw-sst-on? #t) (kw-bsl-on? #f) (kw-std-on? #f) (flow/scheme 24) (pseudo-transient-formulation? #t)");
  CHECK(issue_at(c, "rp/omega/pseudo-relax"));
  CHECK(issue_at(c, "rp/k/pseudo-relax"));
  CHECK_FALSE(issue_at(c, "rp/pressure/relax"));
}

TEST_CASE("Force reports with empty name selections cannot claim resolved zones") {
  auto c = synthetic(config(false, false, true) + R"scm(
    (monitor/report-definitions (((name . "empty")
      (report-definition force-force (surface-names) (force-vector 1 0 0) (scaled? . #f)))))
  )scm");
  REQUIRE(c.reports.size() == 1);
  CHECK_FALSE(c.reports[0].zone_ids);
  CHECK(issue_at(c, "report/empty/zones"));
}

TEST_CASE("Reading a zone selector does not silently accept unknown behavior") {
  for (auto [type, key] : {std::pair{"fluid", "motion-spec"}, {"fluid", "mrf-relative-to-thread"},
                           {"wall", "shear-bc"}, {"wall", "rough-bc"}}) {
    auto c = synthetic(config(false, false, true), "(39 (8 " + std::string(type) +
        " body 1) ((" + key + " . 999)))");
    CHECK(issue_at(c, "zone/8/" + std::string(key)));
  }
}

TEST_CASE("SA case intent is parsed without opening the unverified model allowlist") {
  const auto cfg="(case-config ((rp-3d? . #t) (rp-double? . #t) (rp-unsteady? . #f) (rp-seg? . #t) (rp-lam? . #f) (rp-ke? . #f) (rp-kw? . #f) (rp-sa? . #t))) ";
  auto c=synthetic(std::string(cfg)+"(sa-alternate-prod? #t) (sa-damping? #t) (sa-enhanced-wall-treatment? #t) (sa-cb1 .1355) (sa-mystery 9) (nut/scheme 0)",
    "(39 (4 velocity-inlet inlet 1) ((velocity-spec . 2) (vmag . 1) (sa-spec . 1) (turb-viscosity-ratio . 10) (nut . 3e-5)))");
  CHECK(c.viscous.model==ViscousModel::spalart_allmaras);CHECK(name(*c.viscous.model)=="Spalart-Allmaras");
  CHECK(c.viscous.sa_alternate_production==true);CHECK(c.viscous.sa_damping==true);CHECK(c.viscous.sa_enhanced_wall==true);
  CHECK(c.viscous.coefficients.at("sa-cb1")==doctest::Approx(.1355));
  CHECK(c.numerics.modified_viscosity==Scheme::first_upwind);CHECK_FALSE(c.numerics.k);CHECK_FALSE(c.numerics.omega);
  CHECK(issue_at(c,"viscous-model"));CHECK(issue_at(c,"case-config/rp-sa?"));CHECK(issue_at(c,"rp/sa-mystery"));
  CHECK(boundary(c,4).sa_spec==1);CHECK(boundary(c,4).turbulence.at("nut").constant==doctest::Approx(3e-5));
  auto conflict=synthetic(std::string(cfg).replace(std::string(cfg).find("rp-kw? . #f"),11,"rp-kw? . #t"));
  CHECK_FALSE(conflict.viscous.model);
}
