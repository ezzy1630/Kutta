// CPU reference backend: serial, double-accumulating, deterministic. Every GPU
// backend is tested against it.
#include <algorithm>
#include <cmath>
#include <cstring>
#include <vector>

#include "kernels.hpp"
#include "kutta/runtime/cpu/cpu_kernel_table.hpp"

namespace kutta {
namespace {

class CpuBuffer final : public BufferHandle {
 public:
  explicit CpuBuffer(std::size_t bytes) : data_((bytes + 7) / 8), bytes_(bytes) {}
  std::size_t bytes() const override { return bytes_; }
  void* data() { return data_.data(); }

 private:
  std::vector<std::uint64_t> data_;  // 8-byte aligned storage
  std::size_t bytes_;
};

CpuBuffer& cpu(BufferHandle& h) { return static_cast<CpuBuffer&>(h); }
const float* floats(const Buffer<float>& b) {
  return static_cast<const float*>(cpu(*b.handle()).data());
}

class CpuDevice final : public Device {
 public:
  Backend backend() const override { return Backend::Cpu; }
  std::string name() const override { return "CPU reference (serial)"; }

  double sum(const Buffer<float>& x, std::size_t n) override {
    KUTTA_REQUIRE(n <= x.size(), "sum: n out of range");
    const float* p = floats(x);
    double s = 0;
    for (std::size_t i = 0; i < n; ++i) s += p[i];
    return s;
  }
  double dot(const Buffer<float>& x, const Buffer<float>& y, std::size_t n) override {
    KUTTA_REQUIRE(n <= x.size() && n <= y.size(), "dot: n out of range");
    const float* a = floats(x);
    const float* b = floats(y);
    double s = 0;
    for (std::size_t i = 0; i < n; ++i) s += double(a[i]) * double(b[i]);
    return s;
  }
  double max_abs(const Buffer<float>& x, std::size_t n) override {
    KUTTA_REQUIRE(n <= x.size(), "max_abs: n out of range");
    const float* p = floats(x);
    double m = 0;
    // NaN-propagating, unlike std::max, so a diverged field never reports a finite norm.
    for (std::size_t i = 0; i < n; ++i) {
      const double v = std::fabs(p[i]);
      if (v > m || std::isnan(v)) m = v;
    }
    return m;
  }

  void launch(KernelId id, std::size_t n, std::span<const KernelArg> args) override {
    const auto k = static_cast<std::size_t>(id);
    KUTTA_REQUIRE(k < kKernelCount, "launch: bad kernel id");
    cpu_kernel_table[k](0, n, args.data());
  }
  void synchronize() override {}

 protected:
  void reduce_into(const Buffer<float>& x, const Buffer<float>& y,
      const Buffer<float>& u, const Buffer<float>& v, std::size_t n,
      Buffer<float>& dst, std::size_t offset, Reduction mode) override {
    check_reduction(x, y, u, v, n, dst, offset, mode);
    double a = 0, b = 0;
    for (std::size_t i = 0; i < n; ++i) {
      a += mode == Reduction::Sum ? double(floats(x)[i]) :
          double(floats(x)[i]) * double(floats(y)[i]);
      if (mode == Reduction::Dot2) b += double(floats(u)[i]) * double(floats(v)[i]);
    }
    auto* out = static_cast<float*>(cpu(*dst.handle()).data()) + offset;
    auto store = [](float* p, double value) {
      p[0] = static_cast<float>(value);
      p[1] = std::isfinite(p[0]) ? static_cast<float>(value - double(p[0])) : 0;
    };
    store(out, a);
    if (mode == Reduction::Dot2) store(out + 2, b);
  }
  std::unique_ptr<BufferHandle> alloc_bytes(std::size_t bytes) override {
    return std::make_unique<CpuBuffer>(bytes);
  }
  void upload_bytes(BufferHandle& dst, const void* src, std::size_t bytes) override {
    std::memcpy(cpu(dst).data(), src, bytes);
  }
  void download_bytes(const BufferHandle& src, void* dst, std::size_t bytes) override {
    std::memcpy(dst, cpu(const_cast<BufferHandle&>(src)).data(), bytes);
  }
  void copy_bytes(BufferHandle& dst, const BufferHandle& src, std::size_t bytes) override {
    std::memcpy(cpu(dst).data(), cpu(const_cast<BufferHandle&>(src)).data(), bytes);
  }
  void fill_bytes(BufferHandle& dst, const void* value, std::size_t value_size,
                  std::size_t count) override {
    auto* p = static_cast<unsigned char*>(cpu(dst).data());
    for (std::size_t i = 0; i < count; ++i) std::memcpy(p + i * value_size, value, value_size);
  }
};

}  // namespace

void* cpu_ptr(BufferHandle* h) { return cpu(*h).data(); }

std::unique_ptr<Device> make_cpu_device() { return std::make_unique<CpuDevice>(); }

}  // namespace kutta
