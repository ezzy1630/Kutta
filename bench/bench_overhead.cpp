// Per-iteration overhead on small meshes: dispatch latency, host round trips
// and the V8 realizable k-epsilon channel. Diagnostic only, never a gate.
//
//   bench_overhead latency <cpu|metal|cuda>
//   bench_overhead ke-channel <cpu|metal|cuda> <cells> [counted]
//   bench_overhead parity <cpu|metal|cuda>
//     FNV-1a hashes of every field after fixed iterations of the V8 channel,
//     an SST tetrahedral channel and a laminar box (AMG pressure, automatic
//     pseudo time). Equal hashes across builds mean bitwise-identical results.
#include <chrono>
#include <iomanip>
#include <iostream>
#include <map>

#include "kernels.hpp"
#include "../tests/linalg/instrumented_device.hpp"
#include "../tests/linalg/matrices.hpp"
#include "../validation/cases/ke_channel.hpp"
#include "sst_channel.hpp"

#ifndef KUTTA_BENCH_COMMIT
#define KUTTA_BENCH_COMMIT "unknown"
#endif

namespace {
using Clock = std::chrono::steady_clock;
double seconds_since(Clock::time_point start) {
  return std::chrono::duration<double>(Clock::now() - start).count();
}
kutta::Backend parse_backend(const std::string& name) {
  for (auto b : {kutta::Backend::Cpu, kutta::Backend::Metal, kutta::Backend::Cuda})
    if (name == kutta::backend_name(b)) return b;
  throw kutta::Error("unknown backend: " + name);
}

void latency(kutta::Backend backend) {
  auto dev = kutta::Device::create(backend);
  constexpr std::size_t n = 80;
  auto x = dev->alloc<float>(n), y = dev->alloc<float>(n), out = dev->alloc<float>(4);
  dev->fill(x, 1.0f); dev->fill(y, 0.0f);
  dev->synchronize();
  const auto report = [&](const char* what, int count, double seconds) {
    std::cout << std::setprecision(6) << "OVERHEAD backend=" << kutta::backend_name(backend)
              << " source=" << KUTTA_BENCH_COMMIT << " measure=" << what << " count=" << count
              << " us_each=" << 1e6 * seconds / count << std::endl;
  };
  for (int pass = 0; pass < 2; ++pass) {
    constexpr int dispatches = 20000;
    auto start = Clock::now();
    for (int i = 0; i < dispatches; ++i) kutta::kern::axpy(*dev, n, y, x, 1e-7f);
    dev->synchronize();
    if (pass) report("queued_dispatch", dispatches, seconds_since(start));
    constexpr int reductions = 2000;
    start = Clock::now();
    for (int i = 0; i < reductions; ++i) dev->dot_into(x, y, n, out);
    dev->synchronize();
    if (pass) report("queued_device_reduction", reductions, seconds_since(start));
    {
      auto status = dev->alloc<kutta::Index>(3);
      auto scalars = dev->alloc<float>(24);
      dev->fill(status, 0); dev->fill(scalars, 1.0f);
      auto host_matrix = kutta::test_linalg::poisson(80, 1);
      auto a = kutta::upload(*dev, host_matrix);
      start = Clock::now();
      for (int i = 0; i < reductions; ++i) kutta::kern::solver_bicg_alpha(*dev, 1, scalars, status);
      dev->synchronize();
      if (pass) report("queued_scalar_kernel", reductions, seconds_since(start));
      start = Clock::now();
      for (int i = 0; i < reductions; ++i)
        kutta::kern::csr_spmv_guarded(*dev, n, a.row_ptr, a.col, a.val, x, y, status);
      dev->synchronize();
      if (pass) report("queued_spmv", reductions, seconds_since(start));
      start = Clock::now();
      for (int i = 0; i < reductions; ++i) {
        kutta::kern::csr_spmv_guarded(*dev, n, a.row_ptr, a.col, a.val, x, y, status);
        kutta::kern::solver_bicg_alpha(*dev, 1, scalars, status);
        dev->dot_into(x, y, n, scalars, 2);
        kutta::kern::axpy(*dev, n, y, x, 1e-7f);
      }
      dev->synchronize();
      if (pass) report("queued_mixed_dispatch", 4 * reductions, seconds_since(start));
    }
    constexpr int round_trips = 2000;
    start = Clock::now();
    for (int i = 0; i < round_trips; ++i) { kutta::kern::axpy(*dev, n, y, x, 1e-7f); (void)dev->max_abs(y, n); }
    if (pass) report("dispatch_plus_host_reduction", round_trips, seconds_since(start));
    start = Clock::now();
    for (int i = 0; i < round_trips; ++i) { kutta::kern::axpy(*dev, n, y, x, 1e-7f); dev->synchronize(); }
    if (pass) report("dispatch_plus_synchronize", round_trips, seconds_since(start));
    {
      // Short BiCGStab-Jacobi solves that stay above the float residual floor:
      // per update this is mostly dispatch cost (plus two waits per solve).
      auto host_matrix = kutta::test_linalg::poisson(80, 1, 1, 1, false, false, 0.5);
      auto a = kutta::upload(*dev, host_matrix);
      auto solver = kutta::make_solver(kutta::SolverKind::BicgstabJacobi, *dev, host_matrix);
      auto rhs = dev->alloc<float>(80), solution = dev->alloc<float>(80);
      dev->fill(rhs, 1.0f);
      constexpr int updates = 30, solves = 20;
      int total = 0;
      start = Clock::now();
      for (int k = 0; k < solves; ++k) {
        dev->fill(solution, 0.0f);
        solver->refresh(*dev, a);
        total += solver->solve(*dev, a, solution, rhs, {1e-30, 0, updates, updates}).iterations;
      }
      if (pass) report("bicgstab_jacobi_update", total, seconds_since(start));
    }
    std::vector<float> host(n, 1.0f);
    start = Clock::now();
    for (int i = 0; i < round_trips; ++i) { kutta::kern::axpy(*dev, n, y, x, 1e-7f); dev->upload<float>(x, host); }
    if (pass) report("dispatch_plus_upload", round_trips, seconds_since(start));
  }
}

void ke_channel(kutta::Backend backend, int cells, bool counted) {
  const auto print = [&](const kutta::validation::ke_channel::Result& r, double wall) {
    std::cout << std::setprecision(9) << "KE_CHANNEL backend=" << kutta::backend_name(backend)
              << " source=" << KUTTA_BENCH_COMMIT << " cells=" << cells << " steps=" << r.steps
              << " converged=" << r.converged << " seconds=" << wall
              << " ms_per_step=" << 1e3 * wall / std::max(r.steps, 1) << " tau=" << r.tau
              << " defect=" << r.defect << " change=" << r.change;
  };
  if (!counted) {
    auto dev = kutta::Device::create(backend);
    const auto start = Clock::now();
    const auto r = kutta::validation::ke_channel::run(*dev, cells);
    print(r, seconds_since(start));
    std::cout << std::endl;
    return;
  }
  kutta::test_linalg::InstrumentedDevice dev(backend);
  std::map<std::string, std::size_t> kernels;
  dev.after_launch = [&](kutta::KernelId id, std::size_t, std::span<const kutta::KernelArg>) {
    ++kernels[kutta::kernel_name(id)];
  };
  const auto start = Clock::now();
  const auto r = kutta::validation::ke_channel::run(dev, cells);
  print(r, seconds_since(start));
  const double steps = std::max(r.steps, 1);
  std::cout << " host_waits_per_step=" << dev.host_waits / steps
            << " host_reductions_per_step=" << dev.reductions / steps
            << " downloads_per_step=" << dev.downloads / steps
            << " uploads_per_step=" << dev.uploads / steps
            << " syncs_per_step=" << dev.syncs / steps
            << " launches_per_step=" << dev.launches / steps
            << " device_reductions_per_step=" << dev.device_reductions / steps
            << " allocations_per_step=" << dev.allocations / steps << std::endl;
  for (const auto& [name, calls] : kernels)
    std::cout << "KE_CHANNEL_KERNEL name=" << name << " calls_per_step=" << calls / steps << std::endl;
}
struct Hash {
  std::uint64_t value = 1469598103934665603ull;
  void add(const void* data, std::size_t bytes) {
    const auto* p = static_cast<const unsigned char*>(data);
    for (std::size_t i = 0; i < bytes; ++i) value = (value ^ p[i]) * 1099511628211ull;
  }
  void add(const std::vector<float>& v) { add(v.data(), v.size() * sizeof(float)); }
  void add(double v) { add(&v, sizeof v); }
};

void hash_solver(kutta::Device& dev, const kutta::fv::SteadySolver& solver, bool sst, Hash& h) {
  for (const auto* q : {&solver.velocity().x, &solver.velocity().y, &solver.velocity().z,
                        &solver.pressure(), &solver.mass_flux(), &solver.boundary_pressure()})
    h.add(dev.download(*q));
  if (sst) { h.add(dev.download(solver.kinetic_energy())); h.add(dev.download(solver.specific_dissipation())); }
  const auto& r = solver.report();
  for (const auto* stats : {&r.momentum[0], &r.momentum[1], &r.momentum[2], &r.pressure,
                            &r.turbulence[0], &r.turbulence[1]}) {
    h.add(double(stats->iterations)); h.add(stats->initial_residual); h.add(stats->final_residual);
  }
  for (const auto* v : {&r.u, &r.v, &r.w, &r.continuity, &r.k, &r.omega}) { h.add(v->raw); h.add(v->scaled); }
  h.add(r.pseudo_time_step);
}

void parity(kutta::Backend backend) {
  auto dev = kutta::Device::create(backend);
  const auto print = [&](const char* name, const Hash& h, double seconds) {
    std::cout << "PARITY backend=" << kutta::backend_name(backend) << " source=" << KUTTA_BENCH_COMMIT
              << " case=" << name << " hash=" << std::hex << h.value << std::dec
              << " seconds=" << seconds << std::endl;
  };
  {
    const auto start = Clock::now();
    const auto r = kutta::validation::ke_channel::run(*dev, 20);
    Hash h;
    h.add(double(r.steps)); h.add(r.defect); h.add(r.change); h.add(r.tau);
    for (const auto& p : r.profile) { h.add(p[0]); h.add(p[1]); }
    print("ke-channel-20", h, seconds_since(start));
  }
  {
    const auto start = Clock::now();
    auto mesh = kutta::bench::sst_channel(16, 10, 4);
    const auto bc = kutta::bench::channel_boundaries(mesh);
    auto settings = kutta::bench::channel_settings();
    settings.monitor_interval = 7;
    kutta::fv::SteadySolver solver(*dev, std::move(mesh), bc, settings);
    kutta::bench::initialize_channel(solver);
    solver.iterate(40);
    solver.monitor();
    Hash h;
    hash_solver(*dev, solver, true, h);
    print("sst-tet-channel", h, seconds_since(start));
  }
  {
    const auto start = Clock::now();
    auto mesh = kutta::make_box_hex(12, 12, 6, {0, 0, 0}, {1, 1, 1});
    std::vector<kutta::fv::PatchCondition> bc;
    for (const auto& p : mesh.patches) {
      if (p.name == "zmin" || p.name == "zmax") bc.push_back({p.name, kutta::fv::Symmetry{}});
      else bc.push_back({p.name, kutta::fv::Wall{{p.name == "ymax" ? 1. : 0., 0, 0}}});
    }
    kutta::fv::SolverSettings settings;  // automatic pseudo time step, closed domain
    kutta::fv::SteadySolver solver(*dev, std::move(mesh), bc, settings);
    solver.iterate(30);
    solver.monitor();
    Hash h;
    hash_solver(*dev, solver, false, h);
    print("laminar-lid-box", h, seconds_since(start));
  }
  {
    // Above the Krylov speculation threshold, many reduction groups.
    const auto start = Clock::now();
    auto mesh = kutta::make_box_hex(54, 54, 48, {0, 0, 0}, {1, 1, 1});
    std::vector<kutta::fv::PatchCondition> bc;
    for (const auto& p : mesh.patches)
      bc.push_back({p.name, kutta::fv::Wall{{p.name == "ymax" ? 1. : 0., 0, 0}}});
    kutta::fv::SteadySolver solver(*dev, std::move(mesh), bc, {});
    solver.iterate(3);
    solver.monitor();
    Hash h;
    hash_solver(*dev, solver, false, h);
    print("laminar-lid-box-140k", h, seconds_since(start));
  }
}
}  // namespace

int main(int argc, char** argv) {
  try {
    KUTTA_REQUIRE(argc >= 3, "usage: bench_overhead latency|ke-channel <backend> [cells] [counted]");
    const std::string what = argv[1];
    const auto backend = parse_backend(argv[2]);
    if (what == "latency") latency(backend);
    else if (what == "parity") parity(backend);
    else if (what == "ke-channel") {
      KUTTA_REQUIRE(argc >= 4, "ke-channel needs a cell count");
      ke_channel(backend, std::stoi(argv[3]), argc > 4 && std::string(argv[4]) == "counted");
    } else throw kutta::Error("unknown measurement: " + what);
  } catch (const std::exception& e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
