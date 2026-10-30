// Memory and startup at scale: an SST steady solve on a hex box channel.
// Reports device bytes by array length (per cell, face, matrix entry...),
// host footprint and phase times. Diagnostic only, never a gate.
//
//   bench_memory <cpu|metal|cuda> <nx> <ny> <nz> [iterations]
//   bench_memory <cpu|metal|cuda> --mesh <file.msh|file.cas.h5> [iterations]
// A mesh file keeps its patches: walls, velocity inlets and pressure outlets
// keep their kind, anything else becomes a wall.
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>

#if defined(__APPLE__)
#include <mach/mach.h>
#include <malloc/malloc.h>
#endif
#include <sys/resource.h>
#include <cxxabi.h>
#include <execinfo.h>
#include <cstdlib>

#include "kernels.hpp"
#include "kutta/fv/steady_solver.hpp"
#include "kutta/io/mesh_io.hpp"
#include "../tests/linalg/instrumented_device.hpp"

#ifndef KUTTA_BENCH_COMMIT
#define KUTTA_BENCH_COMMIT "unknown"
#endif

namespace {
using Clock = std::chrono::steady_clock;
double seconds_since(Clock::time_point start) {
  return std::chrono::duration<double>(Clock::now() - start).count();
}
// Process footprint: on Apple silicon this includes shared Metal buffers.
double footprint_gb() {
#if defined(__APPLE__)
  task_vm_info_data_t info{};
  mach_msg_type_number_t count = TASK_VM_INFO_COUNT;
  if (task_info(mach_task_self(), TASK_VM_INFO, reinterpret_cast<task_info_t>(&info), &count) == KERN_SUCCESS)
    return double(info.phys_footprint) / 1e9;
#endif
  return 0;
}
// Highest footprint so far, including transient setup buffers.
double peak_footprint_gb() {
#if defined(__APPLE__)
  task_vm_info_data_t info{};
  mach_msg_type_number_t count = TASK_VM_INFO_COUNT;
  if (task_info(mach_task_self(), TASK_VM_INFO, reinterpret_cast<task_info_t>(&info), &count) == KERN_SUCCESS)
    return double(info.ledger_phys_footprint_peak) / 1e9;
#endif
  return 0;
}
// Host heap in use. Metal buffers are not malloc allocations.
double heap_gb() {
#if defined(__APPLE__)
  malloc_statistics_t stats{};
  malloc_zone_statistics(nullptr, &stats);
  return double(stats.size_in_use) / 1e9;
#else
  return 0;
#endif
}
double peak_rss_gb() {
  rusage usage{};
  getrusage(RUSAGE_SELF, &usage);
#if defined(__APPLE__)
  return double(usage.ru_maxrss) / 1e9;
#else
  return double(usage.ru_maxrss) * 1e3 / 1e9;
#endif
}

// Counts live device allocations by element count, then names the counts
// that match mesh sizes. Wraps InstrumentedDevice so every launch is real.
class SizedDevice final : public kutta::Device {
 public:
  explicit SizedDevice(kutta::Backend b) : inner_(b) {}
  kutta::test_linalg::InstrumentedDevice inner_;
  std::map<std::size_t, std::pair<std::size_t, std::size_t>> live;  // elements -> (arrays, bytes)
  std::map<std::string, std::pair<std::size_t, std::size_t>> sites;  // allocating function -> (arrays, bytes)
  bool attribute = false;
  kutta::Backend backend() const override { return inner_.backend(); }
  std::string name() const override { return inner_.name(); }
  double sum(const kutta::Buffer<float>& x, std::size_t n) override { return inner_.sum(own(x), n); }
  double dot(const kutta::Buffer<float>& x, const kutta::Buffer<float>& y, std::size_t n) override {
    return inner_.dot(own(x), own(y), n);
  }
  double max_abs(const kutta::Buffer<float>& x, std::size_t n) override { return inner_.max_abs(own(x), n); }
  void launch(kutta::KernelId id, std::size_t n, std::span<const kutta::KernelArg> args) override {
    std::vector<kutta::KernelArg> forwarded(args.begin(), args.end());
    for (auto& a : forwarded)
      if (a.kind == kutta::KernelArg::Kind::Buffer) a.buffer = static_cast<Sized*>(a.buffer)->data.handle();
    inner_.launch(id, n, forwarded);
  }
  void synchronize() override { inner_.synchronize(); }

 protected:
  struct Sized final : kutta::BufferHandle {
    kutta::Buffer<float> data;
    SizedDevice* owner;
    std::string site;
    Sized(kutta::Buffer<float> b, SizedDevice* o) : data(std::move(b)), owner(o) {
      auto& [arrays, bytes] = owner->live[data.size()];
      ++arrays; bytes += this->bytes();
      if (owner->attribute) {
        site = caller() + " [" + std::to_string(data.size()) + "]";
        auto& [n, b2] = owner->sites[site];
        ++n; b2 += this->bytes();
      }
    }
    ~Sized() override {
      auto& [arrays, bytes] = owner->live[data.size()];
      --arrays; bytes -= this->bytes();
      if (!site.empty()) { auto& [n, b2] = owner->sites[site]; --n; b2 -= this->bytes(); }
    }
    // First frame outside the Device/Buffer plumbing, demangled, without arguments.
    static std::string caller() {
      void* frames[24];
      const int depth = backtrace(frames, 24);
      char** symbols = backtrace_symbols(frames, depth);
      std::string result = "unknown";
      for (int i = 1; i < depth && symbols; ++i) {
        std::string line = symbols[i];
        const auto mangled_start = line.find(" _Z");
        if (mangled_start == std::string::npos) continue;
        std::string mangled = line.substr(mangled_start + 1);
        mangled = mangled.substr(0, mangled.find(' '));
        int status = 0;
        char* demangled = abi::__cxa_demangle(mangled.c_str(), nullptr, nullptr, &status);
        std::string name = status == 0 && demangled ? demangled : mangled;
        std::free(demangled);
        for (std::size_t at; (at = name.find("(anonymous namespace)")) != std::string::npos;)
          name.replace(at, 21, "anon");
        name = name.substr(0, name.find('('));
        if (name.find("SizedDevice") != std::string::npos || name.find("Device::alloc") != std::string::npos ||
            name.find("Device::upload") != std::string::npos || name.find("make_unique") != std::string::npos ||
            name.find("vector_field") != std::string::npos || name.find("std::") == 0 ||
            name.find("kutta::upload") != std::string::npos)
          continue;
        result = name;
        break;
      }
      std::free(symbols);
      return result;
    }
    std::size_t bytes() const override { return data.size() * sizeof(float); }
  };
  static const kutta::Buffer<float>& own(const kutta::Buffer<float>& b) {
    return static_cast<Sized*>(b.handle())->data;
  }
  void reduce_into(const kutta::Buffer<float>& x, const kutta::Buffer<float>& y, const kutta::Buffer<float>& u,
                   const kutta::Buffer<float>& v, std::size_t n, kutta::Buffer<float>& dst, std::size_t offset,
                   Reduction mode) override {
    auto& out = static_cast<Sized*>(dst.handle())->data;
    if (mode == Reduction::Sum) inner_.sum_into(own(x), n, out, offset);
    else if (mode == Reduction::Dot) inner_.dot_into(own(x), own(y), n, out, offset);
    else inner_.dot2_into(own(x), own(y), own(u), own(v), n, out, offset);
  }
  void host_reductions(std::span<const HostReduction> requests, std::span<double> out) override {
    std::vector<HostReduction> forwarded(requests.begin(), requests.end());
    for (auto& r : forwarded) { r.x = &own(*r.x); r.y = &own(*r.y); }
    const auto values = inner_.reduce_many(forwarded);
    std::copy(values.begin(), values.end(), out.begin());
  }
  std::unique_ptr<kutta::BufferHandle> alloc_bytes(std::size_t bytes) override {
    KUTTA_REQUIRE(bytes % 4 == 0, "memory bench requires 32-bit buffers");
    return std::make_unique<Sized>(inner_.alloc<float>(bytes / 4), this);
  }
  // Forwarded in place: a staging copy here would inflate the footprint
  // being measured. Every buffer holds 32-bit words, and both sides copy bytes.
  void upload_bytes(kutta::BufferHandle& dst, const void* src, std::size_t bytes) override {
    inner_.upload<float>(static_cast<Sized&>(dst).data, {static_cast<const float*>(src), bytes / 4});
  }
  void download_bytes(const kutta::BufferHandle& src, void* dst, std::size_t bytes) override {
    inner_.download<float>(static_cast<const Sized&>(src).data, {static_cast<float*>(dst), bytes / 4});
  }
  void copy_bytes(kutta::BufferHandle& dst, const kutta::BufferHandle& src, std::size_t) override {
    inner_.copy(static_cast<Sized&>(dst).data, static_cast<const Sized&>(src).data);
  }
  void fill_bytes(kutta::BufferHandle& dst, const void* value, std::size_t size, std::size_t) override {
    KUTTA_REQUIRE(size == 4, "memory bench requires 32-bit fills");
    float word;
    std::memcpy(&word, value, 4);
    inner_.fill(static_cast<Sized&>(dst).data, word);
  }
};
}  // namespace

int main(int argc, char** argv) {
  using namespace kutta;
  try {
    KUTTA_REQUIRE(argc >= 4, "usage: bench_memory <cpu|metal|cuda> <nx> <ny> <nz> | --mesh <file> [iterations]");
    const bool from_file = std::string(argv[2]) == "--mesh";
    KUTTA_REQUIRE(from_file || argc >= 5, "bench_memory: box dimensions needed");
    const std::string backend_name_arg = argv[1];
    Backend backend = Backend::Cpu;
    for (auto b : {Backend::Cpu, Backend::Metal, Backend::Cuda})
      if (backend_name_arg == kutta::backend_name(b)) backend = b;
    const int iterations = from_file ? (argc > 4 ? std::stoi(argv[4]) : 3) : (argc > 5 ? std::stoi(argv[5]) : 3);
    std::cout << std::setprecision(6);
    const auto phase = [&](const char* name, double seconds) {
      std::cout << "MEMORY_PHASE source=" << KUTTA_BENCH_COMMIT << " phase=" << name << " seconds=" << seconds
                << " footprint_GB=" << footprint_gb() << " heap_GB=" << heap_gb()
                << " peak_rss_GB=" << peak_rss_gb() << std::endl;
    };
    auto start = Clock::now();
    SizedDevice dev(backend);
    dev.attribute = std::getenv("KUTTA_MEMORY_SITES") != nullptr;
    phase("device", seconds_since(start));
    start = Clock::now();
    auto mesh = from_file ? read_mesh(argv[3]) :
        make_box_hex(std::stoi(argv[2]), std::stoi(argv[3]), std::stoi(argv[4]), {0, 0, 0}, {4, 1, .5});
    std::vector<fv::PatchCondition> bc;
    const fv::TurbulenceValues inflow{.0015, 10};
    if (from_file) for (const auto& p : mesh.patches) {
      if (p.kind == PatchKind::VelocityInlet) bc.push_back({p.name, fv::VelocityInlet{Vec3d{1, 0, 0}, inflow}});
      else if (p.kind == PatchKind::PressureOutlet)
        bc.push_back({p.name, fv::PressureOutlet{0, fv::NormalBackflow{}, fv::BackflowPressure::Static, inflow}});
      else bc.push_back({p.name, fv::Wall{}});
    }
    else for (const auto& p : mesh.patches) {
      if (p.name == "xmin") bc.push_back({p.name, fv::VelocityInlet{Vec3d{1, 0, 0}, inflow}});
      else if (p.name == "xmax")
        bc.push_back({p.name, fv::PressureOutlet{0, fv::NormalBackflow{}, fv::BackflowPressure::Static, inflow}});
      else if (p.name == "ymin" || p.name == "ymax") bc.push_back({p.name, fv::Wall{}});
      else bc.push_back({p.name, fv::Symmetry{}});
    }
    const auto cells = mesh.n_cells, faces = mesh.n_faces(), internal = mesh.n_internal_faces();
    phase("mesh", seconds_since(start));
    fv::SolverSettings settings;
    settings.viscosity = 1e-5; settings.sst = fv::SstSettings{}; settings.pseudo_time_step = .01;
    settings.momentum_relaxation = .3; settings.pressure_relaxation = .2; settings.monitor_interval = 1000;
    start = Clock::now();
    fv::SteadySolver solver(dev, std::move(mesh), bc, settings);
    dev.synchronize();
    phase("solver_setup", seconds_since(start));
    start = Clock::now();
    solver.initialize_fields(std::vector<Vec3d>(cells, {1, 0, 0}), std::vector<double>(cells, 0),
                             std::vector<double>(cells, .0015), std::vector<double>(cells, 10));
    dev.synchronize();
    phase("initialize", seconds_since(start));
    start = Clock::now();
    solver.iterate(1);
    dev.synchronize();
    phase("first_iteration", seconds_since(start));
    start = Clock::now();
    solver.iterate(iterations);
    dev.synchronize();
    const double per_iteration = seconds_since(start) / iterations;
    phase("iterations", per_iteration * iterations);
    const auto nb = std::size_t(faces - internal);
    const auto nnz = std::size_t(cells) + 2 * std::size_t(internal);
    const auto label = [&](std::size_t n) -> std::string {
      if (n == std::size_t(cells)) return "cells";
      if (n == std::size_t(cells) + 1) return "cells+1";
      if (n == std::size_t(faces)) return "faces";
      if (n == std::size_t(internal)) return "internal_faces";
      if (n == nb) return "boundary_faces";
      if (n == nnz) return "matrix_entries";
      if (n == 2 * std::size_t(faces) - nb) return "cell_face_entries";
      return "other";
    };
    std::map<std::string, std::pair<std::size_t, std::size_t>> by_label;
    std::size_t total = 0;
    for (const auto& [elements, use] : dev.live) {
      if (!use.first) continue;
      auto& [arrays, bytes] = by_label[label(elements)];
      arrays += use.first; bytes += use.second; total += use.second;
    }
    for (const auto& [name, use] : by_label)
      std::cout << "MEMORY_ARRAYS source=" << KUTTA_BENCH_COMMIT << " cells=" << cells << " length=" << name
                << " arrays=" << use.first << " bytes_per_cell=" << double(use.second) / cells << std::endl;
    for (const auto& [site, use] : dev.sites)
      if (use.first)
        std::cout << "MEMORY_SITE cells=" << cells << " site=\"" << site << "\" arrays=" << use.first
                  << " bytes_per_cell=" << double(use.second) / cells << std::endl;
    const auto& r = solver.report();
    std::cout << "MEMORY_TOTAL source=" << KUTTA_BENCH_COMMIT << " backend=" << kutta::backend_name(backend)
              << " cells=" << cells << " faces=" << faces << " internal_faces=" << internal
              << " device_bytes_per_cell=" << double(total) / cells
              << " device_peak_bytes_per_cell=" << double(dev.inner_.peak_bytes) / cells
              << " device_GB=" << double(total) / 1e9 << " footprint_GB=" << footprint_gb()
              << " peak_footprint_GB=" << peak_footprint_gb() << " heap_GB=" << heap_gb()
              << " peak_rss_GB=" << peak_rss_gb() << " seconds_per_iteration=" << per_iteration
              << " pressure_iterations=" << r.pressure.iterations
              << " finite=" << std::isfinite(r.pressure.final_residual) << std::endl;
  } catch (const std::exception& e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
