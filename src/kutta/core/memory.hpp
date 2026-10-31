#pragma once

// Process memory readings for setup and solve phases. On Apple silicon the
// footprint includes shared Metal buffers, so host and device are one figure.
// Diagnostic only: nothing here may influence a result.

#include <cstddef>
#include <cstdint>
#include <string_view>

namespace kutta::memory {

struct Usage {
  // Bytes; zero where the platform has no reading.
  std::uint64_t current = 0;  // macOS phys_footprint, Linux RSS, Windows private bytes
  std::uint64_t peak = 0;     // highest `current` since the previous phase mark
  std::uint64_t heap = 0;     // live malloc bytes (macOS only)
};
Usage usage();

// Live device-buffer bytes. Each backend's buffer handle holds one tally.
class DeviceBytes {
 public:
  explicit DeviceBytes(std::size_t bytes);
  ~DeviceBytes();
  DeviceBytes(const DeviceBytes&) = delete;
  DeviceBytes& operator=(const DeviceBytes&) = delete;
  static std::uint64_t live();

 private:
  std::size_t bytes_;
};

// With KUTTA_MEMORY_PHASES set, prints one line per call to stderr:
//   kutta-memory <phase> current=<bytes> peak=<bytes> heap=<bytes> device=<bytes> seconds=<since first mark>
// where peak is the highest footprint reached since the previous mark (exact
// on macOS and Linux; the process lifetime peak on Windows). Otherwise a no-op.
void phase(std::string_view name);

}  // namespace kutta::memory
