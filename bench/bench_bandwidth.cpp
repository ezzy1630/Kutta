#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

#include "kernels.hpp"

#ifndef KUTTA_BENCH_COMMIT
#define KUTTA_BENCH_COMMIT "unknown"
#endif

namespace {
using Clock = std::chrono::steady_clock;
constexpr std::size_t n = 64 * 1024 * 1024;  // 256 MiB per array, beyond the M4 Pro SLC.
constexpr int trials = 7;

template <class F>
double best_seconds(kutta::Device& device, F operation, int repeats = 1) {
  operation();
  device.synchronize();
  double best = std::numeric_limits<double>::infinity();
  for (int trial = 0; trial < trials; ++trial) {
    auto start = Clock::now();
    for (int repeat = 0; repeat < repeats; ++repeat) operation();
    device.synchronize();
    best = std::min(best, std::chrono::duration<double>(Clock::now() - start).count() / repeats);
  }
  return best;
}
std::string json_string(const std::string& value) {
  std::string out = "\"";
  for (char c : value) {
    if (c == '\\' || c == '"') out += '\\';
    if (static_cast<unsigned char>(c) < 32) out += '?';
    else out += c;
  }
  return out + '"';
}
}

int main() {
  try {
    std::cout << "STREAM: 256 MiB/array, 3 arrays, best of " << trials
              << "; decimal GB/s; wall-clock including completion. Commit " << KUTTA_BENCH_COMMIT << '\n';
    std::cout << "backend   copy     scale    add      triad    sum      dot      maxabs   small-sync-us small-batch-us zero-us\n";
    for (auto backend : kutta::Device::available()) {
      auto device = kutta::Device::create(backend);
      auto a = device->alloc<float>(n), b = device->alloc<float>(n), c = device->alloc<float>(n);
      device->fill(a, 1.25f); device->fill(b, .5f); device->fill(c, 0.0f);
      device->synchronize();
      const auto bandwidth = [&](int arrays, auto operation) {
        return double(n * sizeof(float) * arrays) / best_seconds(*device, operation, 8) / 1e9;
      };
      const double copy = bandwidth(2, [&] { kutta::kern::scale_copy(*device, n, c, a, 1.0f); });
      const double scale = bandwidth(2, [&] { kutta::kern::scale_copy(*device, n, c, a, 3.0f); });
      const double add = bandwidth(3, [&] { kutta::kern::axpbyz(*device, n, c, a, b, 1.0f, 1.0f); });
      const double triad = bandwidth(3, [&] { kutta::kern::axpbyz(*device, n, c, a, b, 1.0f, 3.0f); });
      volatile double result = 0;
      const double sum = double(n * sizeof(float)) / best_seconds(*device, [&] { result = device->sum(a, n); }) / 1e9;
      const double dot = double(n * 2 * sizeof(float)) / best_seconds(*device, [&] { result = device->dot(a, b, n); }) / 1e9;
      const double maxabs = double(n * sizeof(float)) / best_seconds(*device, [&] { result = device->max_abs(a, n); }) / 1e9;
      (void)result;
      KUTTA_REQUIRE(device->sum(c, n) == 2.75 * n, "STREAM triad validation failed");
      const double small_sync = best_seconds(*device, [&] {
        kutta::kern::scale_copy(*device, 1, c, a, 1.0f); device->synchronize();
      }, 1000) * 1e6;
      const double small_batch = best_seconds(*device, [&] {
        kutta::kern::scale_copy(*device, 1, c, a, 1.0f);
      }, 1000) * 1e6;
      const double zero = best_seconds(*device, [&] {
        kutta::kern::scale_copy(*device, 0, c, a, 1.0f);
      }, 10000) * 1e6;
      std::cout << std::fixed << std::setprecision(2) << std::left << std::setw(10) << kutta::backend_name(backend);
      for (double value : {copy, scale, add, triad, sum, dot, maxabs, small_sync, small_batch, zero})
        std::cout << std::setw(9) << value;
      std::cout << '\n';
      std::cout << std::setprecision(6)
                << "{\"device\":" << json_string(device->name())
                << ",\"backend\":" << json_string(kutta::backend_name(backend))
                << ",\"commit\":" << json_string(KUTTA_BENCH_COMMIT)
                << ",\"array_bytes\":" << n * sizeof(float) << ",\"trials\":" << trials
                << ",\"dispatches_per_stream_trial\":8"
                << ",\"copy_GB_s\":" << copy << ",\"scale_GB_s\":" << scale
                << ",\"add_GB_s\":" << add << ",\"triad_GB_s\":" << triad
                << ",\"sum_GB_s\":" << sum << ",\"dot_GB_s\":" << dot
                << ",\"max_abs_GB_s\":" << maxabs
                << ",\"small_synchronized_us\":" << small_sync
                << ",\"small_batched_us\":" << small_batch
                << ",\"zero_launch_us\":" << zero << "}\n";
    }
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
