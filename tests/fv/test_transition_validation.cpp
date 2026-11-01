// Intermittency transition model validation (validation/results/transition/
// REPORT.md). Meshes come from `uv run tools/generate_fixtures.py` (never
// committed, found through tests/generated_fixtures.hpp); the restricted data live
// in validation/data and are never committed. Each case skips with a message
// when either is missing, and runs on the first GPU backend present.
//
// Gates are the documented acceptance where Kutta meets it. Where no
// acceptance is documented, or Kutta misses it, the check is a REGRESSION
// GUARD on the measured accuracy, named as such; it is not a pass against the
// reference.
#include <fstream>
#include <map>
#include <sstream>
#include <random>

#include "../backend_helpers.hpp"
#include "../generated_fixtures.hpp"
#include "../../validation/cases/fsae_aero.hpp"

using namespace kutta;
using namespace kutta::validation::fsae_aero;
namespace {
const std::filesystem::path source(KUTTA_SOURCE_DIR);
const std::filesystem::path data = source / "validation/data";
std::unique_ptr<Device> fastest() {
  const auto available = Device::available();
  for (auto b : {Backend::Cuda, Backend::Metal})
    if (std::find(available.begin(), available.end(), b) != available.end()) return Device::create(b);
  return Device::create(Backend::Cpu);
}
bool present(const std::filesystem::path& path, const std::string& how) {
  if (std::filesystem::exists(path)) return true;
  MESSAGE("skipping: missing " << path.string() << " (" << how << ")");
  return false;
}
double number(const std::string& s) {
  std::istringstream in(s);
  in.imbue(std::locale::classic());
  double v = 0;
  in >> v;
  REQUIRE(in);
  return v;
}
std::vector<std::vector<std::string>> csv(const std::filesystem::path& path) {
  std::ifstream in(path);
  REQUIRE(in);
  std::vector<std::vector<std::string>> rows;
  for (std::string line; std::getline(in, line);) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    std::vector<std::string> row;
    std::istringstream cells(line);
    for (std::string cell; std::getline(cells, cell, ',');) {
      if (cell.size() >= 2 && cell.front() == '"' && cell.back() == '"') cell = cell.substr(1, cell.size() - 2);
      row.push_back(cell);
    }
    rows.push_back(row);
  }
  return rows;
}
double interpolate(const std::vector<std::array<double, 2>>& xy, double x) {
  REQUIRE(xy.size() >= 2);
  REQUIRE(x >= xy.front()[0]);
  REQUIRE(x <= xy.back()[0]);
  const auto hi = std::upper_bound(xy.begin(), xy.end(), x, [](double a, const auto& b) { return a < b[0]; });
  if (hi == xy.end()) return xy.back()[1];
  const auto& lo = *(hi - 1);
  return lo[1] + (hi->at(1) - lo[1]) * (x - lo[0]) / (hi->at(0) - lo[0]);
}
std::filesystem::path scratch(const std::string& name) {
  auto path = std::filesystem::temp_directory_path() / ("kutta_" + name + "_" + std::to_string(std::random_device{}()));
  std::filesystem::create_directories(path);
  return path;
}
}  // namespace

// T3A (Coupland, ERCOFTAC case 020): zero-pressure-gradient plate, Tu 3.4% at
// the leading edge. No acceptance for T3A is documented; the experiment bounds
// are the L1 accuracy measured on 2026-10-31 (REPORT.md), rounded up.
// The same solve is also VERIFIED against an independent double-precision
// boundary-layer solution of the documented equations (tracked, our own
// computation: validation/results/transition_verification/REPORT.md).
TEST_CASE("Transition: ERCOFTAC T3A skin friction vs experiment (REGRESSION GUARD) and vs the documented equations (VERIFICATION)" * doctest::test_suite("slow")) {
  const auto summary = data / "ercoftac_t3_flat_plate/raw/t3ay.dat";
  if (!present(summary, "restricted reference data, see its SOURCE.md")) return;
  // Empty, after a skip message, when the generated mesh is missing.
  const auto mesh = generated_fixture("transition/flatplate/t3a_L1.msh");
  if (mesh.empty()) return;
  // Summary rows: tag, x (mm), Re_x, U0, Cf, H, Re_theta, delta, Tu (%).
  std::vector<std::array<double, 2>> measured;  // Re_x, Cf
  std::vector<double> unit;
  std::ifstream in(summary);
  for (std::string line; std::getline(in, line);) {
    std::istringstream row(line);
    std::vector<std::string> t;
    for (std::string s; row >> s;) t.push_back(s);
    if (t.size() != 9 || !t[0].starts_with("1W")) continue;
    measured.push_back({number(t[2]), number(t[4])});
    unit.push_back(number(t[2]) / (number(t[1]) / 1000));
  }
  REQUIRE(measured.size() >= 10);
  std::sort(unit.begin(), unit.end());
  // flatplate.py's Re_unit (median Re_x/x, so U = 1 matches the measured U0)
  // and the inlet turbulence it fits to the measured freestream decay.
  const double re = unit.size() % 2 ? unit[unit.size() / 2] : .5 * (unit[unit.size() / 2 - 1] + unit[unit.size() / 2]);
  Options o;
  o.mesh = mesh;
  o.reynolds = re;
  o.force_walls = {"plate"};
  o.intensity = .046778;
  o.viscosity_ratio = 11.927;
  o.transition = true;
  o.write_cells = true;
  o.output = scratch("t3a");
  // It stops on convergence or its iteration budget, never on wall time.
  auto dev = fastest();
  const auto r = run(*dev, o);
  MESSAGE("T3A " << std::string(backend_name(dev->backend())) << " it=" << r.iterations << " drag=" << r.drag << " residual=" << r.max_residual
                 << " settle_drag=" << r.settle_drag);
  REQUIRE_MESSAGE(r.failure.empty(), r.failure);
  // Drag settles; the continuity residual stalls near 4e-3 on this mesh in
  // laminar, SST and gamma runs alike (REPORT.md), so it is not gated here.
  CHECK(r.settle_drag <= 5e-4);
  // Local edge velocity U_e(x): largest u per cell column below y = 0.1 m;
  // the measured Cf and Re_x use the local freestream speed too.
  std::map<double, double> edge;
  const auto cells = csv(o.output / "cells.csv");
  for (std::size_t i = 1; i < cells.size(); ++i) {
    const double x = number(cells[i][0]), y = number(cells[i][1]), u = number(cells[i][4]);
    if (y >= .1) continue;
    auto [it, inserted] = edge.emplace(x, u);
    if (!inserted) it->second = std::max(it->second, u);
  }
  std::vector<std::array<double, 2>> ue, computed;
  for (const auto& [x, u] : edge) ue.push_back({x, u});
  const auto walls = csv(o.output / "walls.csv");
  for (std::size_t i = 1; i < walls.size(); ++i) {
    if (walls[i][0] != "plate") continue;
    const double x = number(walls[i][1]), u = interpolate(ue, x);
    computed.push_back({x * re * u, -2 * number(walls[i][9]) / (u * u)});
  }
  std::sort(computed.begin(), computed.end());
  std::filesystem::remove_all(o.output);
  // Onset: the Cf minimum between the leading-edge region and full turbulence.
  auto onset = [](const std::vector<std::array<double, 2>>& xy) {
    std::array<double, 2> m{0, 1e30};
    for (const auto& p : xy)
      if (p[0] > 2e4 && p[0] < 4e5 && p[1] < m[1]) m = p;
    return m[0];
  };
  const double onset_error = onset(computed) / onset(measured) - 1;
  double laminar = 0, turbulent = 0;
  for (const auto& [rex, cf] : measured) {
    const double e = interpolate(computed, rex) / cf - 1;
    MESSAGE("Re_x " << rex << " Cf measured " << cf << " Kutta " << interpolate(computed, rex) << " (" << 100 * e << "%)");
    if (rex < onset(measured)) laminar = std::max(laminar, std::abs(e));
    if (rex > 2.6e5) turbulent = std::max(turbulent, std::abs(e));
  }
  MESSAGE("onset Re_x " << onset(computed) << " vs measured " << onset(measured) << " (" << 100 * onset_error << "%)");
  // Measured: onset +11.8%, laminar max 12.1%, turbulent (Re_x > 2.6e5) max 6.9%.
  CHECK(std::abs(onset_error) <= .13);
  CHECK(laminar <= .13);
  CHECK(turbulent <= .075);

  // VERIFICATION: the documented model solved independently
  // (tools/transition/verify/boundary_layer.py, converged to 1.6% in Cf and
  // 0.5% in onset). It predicts the experiment misses above too: Cf 7-15% high
  // before transition. Kutta's L1 mesh is itself within 1-2% of Blasius, so
  // 3% (Cf) and 5% (onset) bound the remaining discretization difference;
  // measured 2026-10-31: 1.7%, 1.3% and +2.6%.
  std::vector<std::array<double, 2>> independent;  // Re_x, Cf gamma
  const auto reference = csv(source / "validation/results/transition_verification/independent_t3a.csv");
  REQUIRE(reference.size() > 100);
  for (std::size_t i = 1; i < reference.size(); ++i) independent.push_back({number(reference[i][1]), number(reference[i][3])});
  const double verified_onset = onset(computed) / onset(independent) - 1;
  double verified_laminar = 0, verified_turbulent = 0;
  for (const auto& [rex, cf] : measured) {
    const double e = interpolate(computed, rex) / interpolate(independent, rex) - 1;
    if (rex < .8 * onset(independent)) verified_laminar = std::max(verified_laminar, std::abs(e));
    if (rex > 2.6e5) verified_turbulent = std::max(verified_turbulent, std::abs(e));
  }
  MESSAGE("vs independent solution: onset " << 100 * verified_onset << "%, laminar max " << 100 * verified_laminar
                                            << "%, turbulent max " << 100 * verified_turbulent << "%");
  CHECK(std::abs(verified_onset) <= .05);
  CHECK(verified_laminar <= .03);
  CHECK(verified_turbulent <= .03);
}

// S1223 at Re 2e5 (UIUC LSATs vol. 1), the documented FSAE acceptance of
// docs/research/verification-cases.md §4: CL within +-0.10 and CD within
// +-25% in the linear range. With transition the drag keeps a +-0.2% wander
// (the laminar separation bubble's), so lift must settle and both
// coefficients are the mean over the last 2000 iterations, as in REPORT.md.
TEST_CASE("Transition: Selig S1223 at Re 2e5, alpha 4.76 deg vs UIUC" * doctest::test_suite("slow")) {
  const auto polars = data / "uiuc_lsat_s1223/s1223_vol1_drag_polars.csv";
  if (!present(polars, "restricted reference data")) return;
  // Empty, after a skip message, when the generated mesh is missing.
  const auto mesh = generated_fixture("fsae_aero/mesh/s1223_a4.76_L2.msh");
  if (mesh.empty()) return;
  double cl = 0, cd = 0;
  for (const auto& row : csv(polars))  // config,re_avg,run_file,alpha_deg,cl,cd
    if (row.size() == 6 && row[0] == "clean" && row[1] == "200400" && row[3] == "4.76") cl = number(row[4]), cd = number(row[5]);
  REQUIRE(cl > 0);
  auto dev = fastest();
  Options o;
  o.mesh = mesh;
  o.reynolds = 2e5;
  o.transition = true;
  // It stops on convergence or its iteration budget, never on wall time.
  const auto r = run(*dev, o);
  REQUIRE_MESSAGE(r.failure.empty(), r.failure);
  double lift = 0, drag = 0, lo = 1e30, hi = -1e30;
  int n = 0;
  for (const auto& s : r.history)
    if (s.iteration > r.iterations - 2000) lift += s.lift, drag += s.drag, lo = std::min(lo, s.drag), hi = std::max(hi, s.drag), ++n;
  REQUIRE(n > 0);
  lift /= n, drag /= n;
  MESSAGE("S1223 gamma " << std::string(backend_name(dev->backend())) << " it=" << r.iterations << " CL=" << lift << " (UIUC " << cl << ") CD=" << drag
                         << " (UIUC " << cd << ") drag band " << lo << ".." << hi << " residual=" << r.max_residual);
  CHECK(r.residuals_met);
  CHECK(r.settle_lift <= 5e-4);
  // Measured 2026-10-31: CL +0.099, CD +16% (fully turbulent SST: +0.075, +23%).
  CHECK(std::abs(lift - cl) <= .10);
  CHECK(std::abs(drag / cd - 1) <= .25);
}
