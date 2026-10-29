#pragma once
// CPU reference backend kernel table, indexed by KernelId (generated).
#include <cstddef>

#include "kutta/runtime/device.hpp"

namespace kutta {

using CpuKernelFn = void (*)(std::size_t begin, std::size_t end, const KernelArg* args);
extern const CpuKernelFn cpu_kernel_table[];

// Host pointer behind a CPU-backend buffer.
void* cpu_ptr(BufferHandle* h);

}  // namespace kutta
