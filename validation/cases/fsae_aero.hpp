#pragma once

// FSAE external-aero validation (docs/research/verification-cases.md, sections
// 4-6): wing sections in and out of ground effect and the Ahmed body. Meshes
// come from tools/fsae_aero/*.py as Fluent .msh files with fixed zone names.
// The kutta CLI only solves Fluent cases, so these runs use the library path,
// as the TMR gates do (validation/cases/tmr.hpp).
#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "kutta/fv/steady_solver.hpp"
#include "kutta/io/mesh_io.hpp"

namespace kutta::validation::fsae_aero {
using namespace fv;

// Nondimensional setup: U = rho = 1 along +x; nu = length/Re.
struct Options {
  std::filesystem::path mesh, output;
  double reynolds = 0, length = 1;   // Re based on `length`
  double reference_area = 1;         // per unit depth in 2D (chord x 1)
  double extrusion = 1;              // 2D meshes are one cell of this depth
  Vec3d lift_direction{0, 1, 0};     // -y reports downforce as positive
  std::vector<std::string> force_walls{"wing"};
  bool moving_ground = false;        // "ground" wall translates at U
  double intensity = .001, viscosity_ratio = 1;
  bool coupled = false;
  // Fluent's intermittency transition model on top of SST, with Fluent's
  // default constants; inflow intermittency 1.
  bool transition = false;
  // Laminar flow (no turbulence model): a mesh check against Blasius.
  bool laminar = false;
  // Fluent hybrid initialization (Laplace potential flow) instead of the
  // uniform freestream; turbulence from the inlet intensity and ratio.
  bool hybrid = false;
  // Convergence-path controls only; they do not change the steady fixed point.
  // Relaxation 0 takes the scheme's Fluent default: SIMPLEC momentum 0.7 and
  // pressure 1.0 (implicit), Coupled 0.5 and 0.5 (pseudo-transient explicit
  // relaxation, as saved by Fluent's Coupled cases such as the V6 elbow).
  double pseudo_time_step = 0, pseudo_time_scale = 1, momentum_relaxation = 0, pressure_relaxation = 0;
  int max_iterations = 6000, check_every = 50;
  bool write_cells = false;  // cells.csv with centroid, velocity, pressure, k, omega, mu_t
  double max_seconds = 3600;
  // Stop once Fluent-scaled residuals are below `residual` AND both force
  // coefficients moved less than `settle` (relative) over `window` iterations.
  double residual = 1e-3, settle = 5e-4;
  int window = 500;
};
struct Sample {
  int iteration = 0;
  double seconds = 0, lift = 0, drag = 0, drag_pressure = 0;
  IterationReport report;
};
struct Result {
  Index cells = 0;
  int iterations = 0;
  double setup_seconds = 0, solve_seconds = 0, lift = 0, drag = 0, drag_pressure = 0;
  double settle_lift = 0, settle_drag = 0, max_residual = 0;
  double wall_yplus_max = 0, wall_yplus_mean = 0;
  bool residuals_met = false, settled = false, stopped_for_budget = false;
  std::string failure;
  std::vector<Sample> history;
};

inline double max_residual(const IterationReport& r, bool three_d, bool transition = false) {
  double m = std::max({r.u.fluent_scaled, r.v.fluent_scaled, r.continuity.fluent_scaled, r.k.fluent_scaled, r.omega.fluent_scaled});
  // Laminar runs report zero turbulence residuals.
  if (transition) m = std::max(m, r.intermittency.fluent_scaled);
  return three_d ? std::max(m, r.w.fluent_scaled) : m;
}

inline std::vector<PatchCondition> boundaries(const PolyMesh& mesh, const Options& o) {
  const TurbulenceCondition inflow = IntensityViscosityRatio{o.intensity, o.viscosity_ratio};
  std::vector<PatchCondition> bc;
  for (const auto& p : mesh.patches) {
    if (p.kind == PatchKind::Wall) {
      Wall wall;
      if (p.name == "ground" && o.moving_ground) wall.velocity = Vec3d{1, 0, 0};
      bc.push_back({p.name, wall});
    } else if (p.kind == PatchKind::VelocityInlet) bc.push_back({p.name, VelocityInlet{Vec3d{1, 0, 0}, inflow}});
    else if (p.kind == PatchKind::PressureOutlet) bc.push_back({p.name, PressureOutlet{0, NormalBackflow{}, BackflowPressure::Static, inflow}});
    else if (p.kind == PatchKind::Symmetry) bc.push_back({p.name, Symmetry{}});
    else throw Error("fsae_aero: unmapped boundary " + p.name + " (" + p.source_type + ")");
  }
  return bc;
}

inline Sample sample(SteadySolver& solver, const Options& o, int iteration, double seconds) {
  const auto f = solver.forces(o.force_walls);
  const ReferenceValues ref{1, 1, o.reference_area, 0};
  Sample s;
  s.iteration = iteration;
  s.seconds = seconds;
  s.lift = f.coefficient(o.lift_direction, ref);
  s.drag = f.coefficient({1, 0, 0}, ref);
  ForceReport pressure_only = f;
  pressure_only.viscous = {};
  s.drag_pressure = pressure_only.coefficient({1, 0, 0}, ref);
  s.report = solver.report();
  return s;
}

// Wall distributions: one row per wall face in the original patch order.
inline void write_walls(Device& dev, SteadySolver& solver, const Options& o, std::ostream& out, Result& r) {
  const auto pressure = dev.download(solver.boundary_pressure());
  const auto& shear = solver.wall_shear();
  const auto tx = dev.download(shear.x), ty = dev.download(shear.y), tz = dev.download(shear.z);
  const auto yplus = o.laminar ? std::vector<float>(tx.size(), 0.0f) : dev.download(solver.turbulence().wall_y_plus);
  const auto& mesh = solver.mesh();
  const auto& g = solver.geometry();
  out << "patch,x,y,z,nx,ny,nz,area,cp,tau_x,tau_y,tau_z,yplus\n" << std::setprecision(9);
  double sum = 0;
  Index count = 0;
  for (const auto& p : mesh.patches) {
    if (p.kind != PatchKind::Wall) continue;
    const bool force = std::find(o.force_walls.begin(), o.force_walls.end(), p.name) != o.force_walls.end();
    for (Index f = p.start; f < p.start + p.size; ++f) {
      const auto b = f - mesh.n_internal_faces();
      const auto a = g.face_area[f];
      const double area = std::sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
      const auto c = g.face_centroid[f];
      out << p.name << ',' << c.x << ',' << c.y << ',' << c.z << ',' << a.x / area << ',' << a.y / area << ',' << a.z / area << ','
          << area << ',' << 2 * pressure[b] << ',' << tx[b] << ',' << ty[b] << ',' << tz[b] << ',' << yplus[b] << '\n';
      if (force) {
        r.wall_yplus_max = std::max(r.wall_yplus_max, double(yplus[b]));
        sum += yplus[b];
        ++count;
      }
    }
  }
  r.wall_yplus_mean = count ? sum / count : 0;
}

inline double relative_change(const std::vector<Sample>& h, int window, double Sample::*value) {
  const auto& last = h.back();
  double change = 0;
  for (auto it = h.rbegin(); it != h.rend() && last.iteration - it->iteration <= window; ++it)
    change = std::max(change, std::abs((*it).*value - last.*value) / std::max(std::abs(last.*value), 1e-12));
  return h.front().iteration > last.iteration - window ? std::numeric_limits<double>::infinity() : change;
}

inline Result run(Device& dev, const Options& o) {
  using Clock = std::chrono::steady_clock;
  auto start = Clock::now();
  KUTTA_REQUIRE(std::filesystem::exists(o.mesh), "fsae_aero mesh missing (generate with tools/fsae_aero): " + o.mesh.string());
  KUTTA_REQUIRE(o.reynolds > 0, "fsae_aero: Reynolds number required");
  auto mesh = read_fluent_msh(o.mesh, o.extrusion);
  const bool three_d = std::none_of(mesh.patches.begin(), mesh.patches.end(), [](const Patch& p) { return p.name == "front"; });
  auto bc = boundaries(mesh, o);
  SolverSettings s;
  s.viscosity = o.length / o.reynolds;
  KUTTA_REQUIRE(!(o.laminar && o.transition), "fsae_aero: transition requires the turbulence model");
  if (!o.laminar) s.sst = SstSettings{};
  if (o.transition) s.sst->transition = IntermittencySettings{};
  s.monitor_interval = o.check_every;
  if (o.coupled) s.coupling = PressureVelocityCoupling::Coupled;
  s.pseudo_time_step = o.pseudo_time_step;
  s.pseudo_time_scale = o.pseudo_time_scale;
  s.momentum_relaxation = o.momentum_relaxation > 0 ? o.momentum_relaxation : o.coupled ? .5 : .7;
  s.pressure_relaxation = o.pressure_relaxation > 0 ? o.pressure_relaxation : o.coupled ? .5 : 1;
  SteadySolver solver(dev, std::move(mesh), std::move(bc), s);
  const auto n = solver.mesh().n_cells;
  const double k = 1.5 * o.intensity * o.intensity, omega = k / (s.viscosity * o.viscosity_ratio);
  std::vector<Vec3d> u(n, Vec3d{1, 0, 0});
  std::vector<double> p(n, 0), kk(n, k), w(n, omega);
  if (o.hybrid) {
    HybridInitialization h;
    h.turbulence_intensity = o.intensity;
    h.turbulence_viscosity_ratio = o.viscosity_ratio;
    solver.initialize_hybrid(h);
  } else {
    if (o.laminar) solver.initialize_fields(u, p);
    else solver.initialize_fields(u, p, kk, w);
  }
  dev.synchronize();
  Result r;
  r.cells = n;
  r.setup_seconds = std::chrono::duration<double>(Clock::now() - start).count();
  start = Clock::now();
  std::ofstream history;
  if (!o.output.empty()) {
    std::filesystem::create_directories(o.output);
    history.open(o.output / "history.csv");
    history << "iteration,seconds,lift,drag,drag_pressure,u,v,w,continuity,k,omega," << (o.transition ? "intermit," : "")
            << "max_velocity,continuity_raw\n" << std::setprecision(10);
  }
  while (r.iterations < o.max_iterations) {
    try {
      solver.iterate(std::min(o.check_every, o.max_iterations - r.iterations));
      solver.monitor();
      r.history.push_back(sample(solver, o, solver.report().iteration, std::chrono::duration<double>(Clock::now() - start).count()));
    } catch (const Error& e) {
      r.failure = e.what();
      break;
    }
    const auto& h = r.history.back();
    r.iterations = h.iteration;
    r.lift = h.lift;
    r.drag = h.drag;
    r.drag_pressure = h.drag_pressure;
    r.max_residual = max_residual(h.report, three_d, o.transition);
    r.settle_lift = relative_change(r.history, o.window, &Sample::lift);
    r.settle_drag = relative_change(r.history, o.window, &Sample::drag);
    r.residuals_met = r.max_residual <= o.residual;
    r.settled = r.settle_lift <= o.settle && r.settle_drag <= o.settle;
    if (history) {
      const auto& q = h.report;
      history << h.iteration << ',' << h.seconds << ',' << h.lift << ',' << h.drag << ',' << h.drag_pressure << ',' << q.u.fluent_scaled << ','
              << q.v.fluent_scaled << ',' << q.w.fluent_scaled << ',' << q.continuity.fluent_scaled << ',' << q.k.fluent_scaled << ','
              << q.omega.fluent_scaled << ',';
      if (o.transition) history << q.intermittency.fluent_scaled << ',';
      // Raw continuity compares schemes: the Fluent scale is the first five
      // iterations' imbalance, which Coupled keeps far smaller than SIMPLEC.
      history << q.maximum_velocity << ',' << q.continuity.raw << '\n' << std::flush;
    }
    std::cout << std::setprecision(8) << "fsae_aero it=" << r.iterations << " lift=" << r.lift << " drag=" << r.drag
              << " residual=" << r.max_residual << " settle=" << std::max(r.settle_lift, r.settle_drag) << std::endl;
    if (!std::isfinite(r.lift) || !std::isfinite(r.drag) || !std::isfinite(r.max_residual)) {
      r.failure = "non-finite forces or residuals";
      break;
    }
    if (r.residuals_met && r.settled) break;
    if (std::chrono::duration<double>(Clock::now() - start).count() >= o.max_seconds) {
      r.stopped_for_budget = true;
      break;
    }
  }
  dev.synchronize();
  if (o.write_cells && !o.output.empty()) {
    const auto ux = dev.download(solver.velocity().x), uy = dev.download(solver.velocity().y), uz = dev.download(solver.velocity().z);
    // Laminar runs have no turbulence fields; their columns are zero.
    const auto turbulent = [&](auto field) { return o.laminar ? std::vector<float>(std::size_t(n), 0.0f) : dev.download(field()); };
    const auto pc = dev.download(solver.pressure());
    const auto kc = turbulent([&]() -> const Buffer<float>& { return solver.kinetic_energy(); });
    const auto wc = turbulent([&]() -> const Buffer<float>& { return solver.specific_dissipation(); });
    const auto mt = turbulent([&]() -> const Buffer<float>& { return solver.turbulence().mu_t; });
    const auto gamma = o.transition ? dev.download(solver.intermittency()) : std::vector<float>{};
    std::ofstream cells(o.output / "cells.csv");
    cells << "x,y,z,volume,u,v,w,p,k,omega,mu_t" << (o.transition ? ",intermittency" : "") << "\n" << std::setprecision(9);
    for (Index i = 0; i < n; ++i) {
      const auto c = solver.geometry().cell_centroid[i];
      cells << c.x << ',' << c.y << ',' << c.z << ',' << solver.geometry().cell_volume[i] << ',' << ux[i] << ',' << uy[i] << ',' << uz[i] << ','
            << pc[i] << ',' << kc[i] << ',' << wc[i] << ',' << mt[i];
      if (o.transition) cells << ',' << gamma[i];
      cells << '\n';
    }
  }
  r.solve_seconds = std::chrono::duration<double>(Clock::now() - start).count();
  if (r.failure.empty()) {
    if (!o.output.empty()) {
      std::ofstream walls(o.output / "walls.csv");
      write_walls(dev, solver, o, walls, r);
    } else {
      std::ostringstream discard;
      write_walls(dev, solver, o, discard, r);
    }
  }
  return r;
}

inline void print(std::ostream& out, const Result& r) {
  // JSON has no infinity; unsettled histories print null.
  const auto num = [&](double v) -> std::ostream& { return std::isfinite(v) ? out << v : out << "null"; };
  out << std::setprecision(10) << "{\"cells\":" << r.cells << ",\"iterations\":" << r.iterations << ",\"setup_s\":" << r.setup_seconds
      << ",\"solve_s\":" << r.solve_seconds << ",\"lift\":";
  num(r.lift) << ",\"drag\":";
  num(r.drag) << ",\"drag_pressure\":";
  num(r.drag_pressure) << ",\"max_residual\":";
  num(r.max_residual) << ",\"settle_lift\":";
  num(r.settle_lift) << ",\"settle_drag\":";
  num(r.settle_drag) << ",\"residuals_met\":" << (r.residuals_met ? "true" : "false") << ",\"settled\":" << (r.settled ? "true" : "false")
      << ",\"budget\":" << (r.stopped_for_budget ? "true" : "false") << ",\"yplus_max\":" << r.wall_yplus_max
      << ",\"yplus_mean\":" << r.wall_yplus_mean << ",\"failure\":\"";
  for (char c : r.failure) out << (c == '"' || c == '\\' ? ' ' : c == '\n' ? ';' : c);
  out << "\"}";
}
}  // namespace kutta::validation::fsae_aero
