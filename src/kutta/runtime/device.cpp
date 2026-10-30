#include "kutta/runtime/device.hpp"
#if defined(_WIN32)
#define NOMINMAX
#include <windows.h>
#elif defined(__APPLE__)
#include <sys/sysctl.h>
#else
#include <unistd.h>
#endif

namespace kutta {

std::unique_ptr<Device> make_cpu_device();
#if KUTTA_HAS_METAL
std::unique_ptr<Device> make_metal_device();  // throws if no Metal device
std::vector<DeviceInfo> enumerate_metal_devices();
#endif
#if KUTTA_HAS_CUDA
std::unique_ptr<Device> make_cuda_device();  // throws if no CUDA device
std::vector<DeviceInfo> enumerate_cuda_devices();
#endif

const char* backend_name(Backend b) {
  switch (b) {
    case Backend::Cpu: return "cpu";
    case Backend::Metal: return "metal";
    case Backend::Cuda: return "cuda";
  }
  return "?";
}

std::unique_ptr<Device> Device::create(Backend b) {
  switch (b) {
    case Backend::Cpu: return make_cpu_device();
    case Backend::Metal:
#if KUTTA_HAS_METAL
      return make_metal_device();
#else
      throw Error("Kutta was built without the Metal backend");
#endif
    case Backend::Cuda:
#if KUTTA_HAS_CUDA
      return make_cuda_device();
#else
      throw Error("Kutta was built without the CUDA backend");
#endif
  }
  throw Error("unknown backend");
}

std::vector<Backend> Device::available() {
  std::vector<Backend> out{Backend::Cpu};
  for (Backend b : {Backend::Metal, Backend::Cuda}) {
    try {
      (void)create(b);
      out.push_back(b);
    } catch (const Error&) {
    }
  }
  return out;
}

std::vector<DeviceInfo> Device::enumerate() {
  std::vector<DeviceInfo> out;
#if KUTTA_HAS_CUDA
  out = enumerate_cuda_devices();
#endif
#if KUTTA_HAS_METAL
  auto metal = enumerate_metal_devices();
  for (auto& d : metal) {
    if (!out.empty()) d.is_default = false;
  }
  out.insert(out.end(), metal.begin(), metal.end());
#endif
  std::uint64_t memory = 0;
#if defined(_WIN32)
  MEMORYSTATUSEX state{}; state.dwLength = sizeof(state);
  if (GlobalMemoryStatusEx(&state)) memory = state.ullTotalPhys;
#elif defined(__APPLE__)
  std::size_t size = sizeof(memory);
  if (sysctlbyname("hw.memsize", &memory, &size, nullptr, 0) != 0) memory = 0;
#else
  const auto pages = sysconf(_SC_PHYS_PAGES), page_size = sysconf(_SC_PAGESIZE);
  if (pages > 0 && page_size > 0) memory = std::uint64_t(pages) * std::uint64_t(page_size);
#endif
  out.push_back({Backend::Cpu, "CPU reference (serial)", memory, 0, out.empty()});
  return out;
}

}  // namespace kutta
