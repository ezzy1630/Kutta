#pragma once

namespace kutta::fv {
// Recompute trades operator throughput for memory on capacity-limited devices.
// Both paths preserve the same float coefficient rounding boundaries.
enum class CoupledEdgeStorage { Cached, Recompute };
}  // namespace kutta::fv
