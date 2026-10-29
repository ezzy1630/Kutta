#pragma once

#include <algorithm>
#include <cstring>
#include <functional>

#include "kutta/runtime/device.hpp"

namespace kutta::test_linalg {
// Run real kernels/reductions on the selected backend while counting public
// Device costs. The wrapper never emulates numerical work.
class InstrumentedDevice final : public Device {
  struct Allocation final : BufferHandle {
    Buffer<float> data;
    std::size_t& live;
    Allocation(Buffer<float> b, std::size_t& count) : data(std::move(b)), live(count) {
      live += bytes();
    }
    ~Allocation() override { live -= bytes(); }
    std::size_t bytes() const override { return data.size() * sizeof(float); }
  };
  static Allocation& allocation(BufferHandle& h) { return static_cast<Allocation&>(h); }
  static const Allocation& allocation(const BufferHandle& h) {
    return static_cast<const Allocation&>(h);
  }
 public:
  explicit InstrumentedDevice(Backend backend) : inner_(Device::create(backend)) {}
  Backend backend() const override { return inner_->backend(); }
  std::string name() const override { return inner_->name(); }
  std::size_t live_bytes = 0, peak_bytes = 0, reductions = 0, downloads = 0, device_reductions = 0;
  std::function<void(KernelId, std::size_t, std::span<const KernelArg>)> after_launch;
  double sum(const Buffer<float>& x, std::size_t n) override {
    ++reductions; return inner_->sum(allocation(*x.handle()).data, n);
  }
  double dot(const Buffer<float>& x, const Buffer<float>& y, std::size_t n) override {
    ++reductions;
    return inner_->dot(allocation(*x.handle()).data, allocation(*y.handle()).data, n);
  }
  double max_abs(const Buffer<float>& x, std::size_t n) override {
    ++reductions; return inner_->max_abs(allocation(*x.handle()).data, n);
  }
  void launch(KernelId id, std::size_t n, std::span<const KernelArg> args) override {
    std::vector<KernelArg> forwarded(args.begin(), args.end());
    for (auto& arg : forwarded)
      if (arg.kind == KernelArg::Kind::Buffer) arg.buffer = allocation(*arg.buffer).data.handle();
    inner_->launch(id, n, forwarded);
    if (after_launch) after_launch(id, n, args);
  }
  template<class T> std::vector<T> read(BufferHandle& h) {
    auto words = inner_->download(allocation(h).data);
    std::vector<T> result(words.size());
    if (!words.empty()) std::memcpy(result.data(), words.data(), words.size() * sizeof(float));
    return result;
  }
  void synchronize() override { inner_->synchronize(); }
 protected:
  void reduce_into(const Buffer<float>& x, const Buffer<float>& y,
      const Buffer<float>& u, const Buffer<float>& v, std::size_t n,
      Buffer<float>& dst, std::size_t offset, Reduction mode) override {
    check_reduction(x, y, u, v, n, dst, offset, mode);
    ++device_reductions;
    auto& out = allocation(*dst.handle()).data;
    const auto& a = allocation(*x.handle()).data;
    const auto& b = allocation(*y.handle()).data;
    if (mode == Reduction::Sum) inner_->sum_into(a, n, out, offset);
    else if (mode == Reduction::Dot) inner_->dot_into(a, b, n, out, offset);
    else inner_->dot2_into(a, b, allocation(*u.handle()).data,
                           allocation(*v.handle()).data, n, out, offset);
  }
  std::unique_ptr<BufferHandle> alloc_bytes(std::size_t bytes) override {
    KUTTA_REQUIRE(bytes % sizeof(float) == 0, "instrumentation requires 32-bit buffers");
    auto result = std::make_unique<Allocation>(inner_->alloc<float>(bytes / sizeof(float)), live_bytes);
    peak_bytes = std::max(peak_bytes, live_bytes);
    return result;
  }
  void upload_bytes(BufferHandle& dst, const void* src, std::size_t bytes) override {
    std::vector<float> words(bytes / sizeof(float));
    if (bytes) std::memcpy(words.data(), src, bytes);
    inner_->upload<float>(allocation(dst).data, words);
  }
  void download_bytes(const BufferHandle& src, void* dst, std::size_t bytes) override {
    ++downloads;
    auto words = inner_->download(allocation(src).data);
    if (bytes) std::memcpy(dst, words.data(), bytes);
  }
  void copy_bytes(BufferHandle& dst, const BufferHandle& src, std::size_t) override {
    inner_->copy(allocation(dst).data, allocation(src).data);
  }
  void fill_bytes(BufferHandle& dst, const void* value, std::size_t size, std::size_t) override {
    KUTTA_REQUIRE(size == sizeof(float), "instrumentation requires 32-bit fills");
    float word; std::memcpy(&word, value, size);
    inner_->fill(allocation(dst).data, word);
  }
 private:
  std::unique_ptr<Device> inner_;
};
}  // namespace kutta::test_linalg
