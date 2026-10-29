#include "kutta/runtime/device.hpp"

namespace kutta {

std::unique_ptr<Device> make_cpu_device();
#if KUTTA_HAS_METAL
std::unique_ptr<Device> make_metal_device();  // throws if no Metal device
#endif
#if KUTTA_HAS_CUDA
std::unique_ptr<Device> make_cuda_device();  // throws if no CUDA device
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

}  // namespace kutta
