#include "kutta/core/memory.hpp"

#include <atomic>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <mutex>
#include <string>

#if defined(__APPLE__)
#include <libproc.h>
#include <malloc/malloc.h>
#include <sys/resource.h>
#include <unistd.h>
// Exported by libSystem; declared only in the private headers.
extern "C" int proc_rlimit_control(pid_t pid, int flavor, void* arg);
#elif defined(_WIN32)
#define NOMINMAX
#include <windows.h>
#include <psapi.h>
#endif

namespace kutta::memory {
namespace {
std::atomic<std::uint64_t> device_live{0};
#if defined(__linux__)
std::uint64_t status_kib(const std::string& key) {
  std::ifstream in("/proc/self/status");
  for (std::string line; std::getline(in, line);)
    if (line.rfind(key, 0) == 0) return std::strtoull(line.c_str() + key.size(), nullptr, 10) * 1024;
  return 0;
}
#endif
void reset_peak() {
#if defined(__APPLE__)
  proc_rlimit_control(getpid(), RLIMIT_FOOTPRINT_INTERVAL, reinterpret_cast<void*>(FOOTPRINT_INTERVAL_RESET));
#elif defined(__linux__)
  std::ofstream("/proc/self/clear_refs") << "5";
#endif
}
}  // namespace

Usage usage() {
#if defined(__APPLE__)
  rusage_info_v4 info{};
  malloc_statistics_t heap{};
  malloc_zone_statistics(nullptr, &heap);
  if (proc_pid_rusage(getpid(), RUSAGE_INFO_V4, reinterpret_cast<rusage_info_t*>(&info)) == 0)
    return {info.ri_phys_footprint, info.ri_interval_max_phys_footprint, heap.size_in_use};
#elif defined(__linux__)
  return {status_kib("VmRSS:"), status_kib("VmHWM:")};
#elif defined(_WIN32)
  PROCESS_MEMORY_COUNTERS_EX counters{};
  if (GetProcessMemoryInfo(GetCurrentProcess(), reinterpret_cast<PROCESS_MEMORY_COUNTERS*>(&counters), sizeof counters))
    return {counters.PrivateUsage, counters.PeakPagefileUsage};
#endif
  return {};
}

DeviceBytes::DeviceBytes(std::size_t bytes) : bytes_(bytes) { device_live += bytes; }
DeviceBytes::~DeviceBytes() { device_live -= bytes_; }
std::uint64_t DeviceBytes::live() { return device_live; }

void phase(std::string_view name) {
  static const bool enabled = std::getenv("KUTTA_MEMORY_PHASES") != nullptr;
  if (!enabled) return;
  static std::mutex lock;
  static const auto start = std::chrono::steady_clock::now();
  const std::lock_guard guard(lock);
  const auto u = usage();
  const auto seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
  std::fprintf(stderr, "kutta-memory %.*s current=%llu peak=%llu heap=%llu device=%llu seconds=%.3f\n",
               int(name.size()), name.data(), static_cast<unsigned long long>(u.current),
               static_cast<unsigned long long>(u.peak), static_cast<unsigned long long>(u.heap),
               static_cast<unsigned long long>(DeviceBytes::live()), seconds);
  std::fflush(stderr);
  reset_peak();
}

}  // namespace kutta::memory
