#import <Foundation/Foundation.h>
#import <Metal/Metal.h>

#include <algorithm>
#include <climits>
#include <cmath>
#include <cstring>
#include <vector>

#include "kernels.hpp"
#include "kutta/runtime/metal/metal_kernel_table.hpp"

namespace kutta {
namespace {

// TwoSum keeps the rounding residual at every addition, including across SIMD
// lanes and threadgroups. Dot products also retain the FMA product residual.
// With fast math disabled these float pairs have roughly 48 bits of precision;
// the host finishes the bounded set of (hi, lo) partials in double. Range is
// still float's: partial sums beyond ~3e38 become Inf, and products below ~1e-38
// flush to zero, so callers keep solver vectors O(1)-scaled. The max is
// NaN-propagating so a diverged field can never report a finite norm.
constexpr const char* reduction_source = R"MSL(
#include <metal_stdlib>
using namespace metal;
float2 pair_add(float2 a, float2 b) {
  float s = a.x + b.x;
  if (!isfinite(s)) return float2(s, 0.0f);  // keep Inf as Inf instead of Inf - Inf = NaN
  float v = s - a.x;
  float e = (a.x - (s - v)) + (b.x - v);
  e = e + a.y + b.y;
  float h = s + e;
  return float2(h, e - (h - s));
}
float nan_max(float a, float b) { return (b > a || isnan(b)) ? b : a; }
kernel void runtime_reduce(const device float* x [[buffer(0)]],
                           const device float* y [[buffer(1)]],
                           device float2* partial [[buffer(2)]],
                           constant uint& n [[buffer(3)]],
                           constant uint& mode [[buffer(4)]],
                           const device float* u [[buffer(5)]],
                           const device float* v [[buffer(6)]],
                           uint gid [[thread_position_in_grid]],
                           uint tid [[thread_index_in_threadgroup]],
                           uint group [[threadgroup_position_in_grid]],
                           uint width [[threads_per_simdgroup]],
                           uint lane [[thread_index_in_simdgroup]],
                           uint simd [[simdgroup_index_in_threadgroup]],
                           uint threads [[threads_per_threadgroup]],
                           uint groups [[threadgroups_per_grid]]) {
  threadgroup float2 scratch[32];
  threadgroup float2 second[32];
  float2 a = float2(0.0f), c = float2(0.0f);
  for (uint i = gid; i < n; i += threads * groups) {
    if (mode == 2) a.x = nan_max(a.x, abs(x[i]));
    else {
      float h = mode == 1 || mode == 3 ? x[i] * y[i] : x[i];
      float l = mode == 1 || mode == 3 ? fma(x[i], y[i], -h) : 0.0f;
      a = pair_add(a, float2(h, l));
      if (mode == 3) {
        float h2 = u[i] * v[i];
        c = pair_add(c, float2(h2, fma(u[i], v[i], -h2)));
      }
    }
  }
  for (uint d = width / 2; d > 0; d /= 2) {
    float2 b = simd_shuffle_down(a, d);
    float2 e = simd_shuffle_down(c, d);
    if (mode == 3 && lane + d < width) c = pair_add(c, e);
    if (lane + d < width) a = mode == 2 ? float2(nan_max(a.x, b.x), 0) : pair_add(a, b);
  }
  if (lane == 0) { scratch[simd] = a; second[simd] = c; }
  threadgroup_barrier(mem_flags::mem_threadgroup);
  if (simd == 0) {
    a = lane < threads / width ? scratch[lane] : float2(0.0f);
    c = lane < threads / width ? second[lane] : float2(0.0f);
    for (uint d = width / 2; d > 0; d /= 2) {
      float2 b = simd_shuffle_down(a, d);
      float2 e = simd_shuffle_down(c, d);
      if (mode == 3 && lane + d < width) c = pair_add(c, e);
      if (lane + d < width) a = mode == 2 ? float2(nan_max(a.x, b.x), 0) : pair_add(a, b);
    }
    if (lane == 0) {
      partial[group] = a;
      if (mode == 3) partial[1024 + group] = c;
    }
  }
}
// A second fixed tree finishes partials on the device. A serial 1024-pair
// finish would replace the host wait with a long single-lane dependency chain.
kernel void runtime_reduce_into(const device float2* partial [[buffer(0)]],
                                device float* out [[buffer(1)]],
                                constant uint& groups [[buffer(2)]],
                                constant uint& offset [[buffer(3)]],
                                constant uint& dual [[buffer(4)]],
                                uint tid [[thread_index_in_threadgroup]],
                                uint width [[threads_per_simdgroup]],
                                uint lane [[thread_index_in_simdgroup]],
                                uint simd [[simdgroup_index_in_threadgroup]],
                                uint threads [[threads_per_threadgroup]]) {
  threadgroup float2 first[32], second[32];
  float2 a = float2(0), b = float2(0);
  for (uint j = tid; j < groups; j += threads) {
    a = pair_add(a, partial[j]);
    if (dual) b = pair_add(b, partial[1024 + j]);
  }
  for (uint d = width / 2; d > 0; d /= 2) {
    float2 c = simd_shuffle_down(a, d), e = simd_shuffle_down(b, d);
    if (lane + d < width) { a = pair_add(a, c); if (dual) b = pair_add(b, e); }
  }
  if (lane == 0) { first[simd] = a; second[simd] = b; }
  threadgroup_barrier(mem_flags::mem_threadgroup);
  if (simd == 0) {
    a = lane < threads / width ? first[lane] : float2(0);
    b = lane < threads / width ? second[lane] : float2(0);
    for (uint d = width / 2; d > 0; d /= 2) {
      float2 c = simd_shuffle_down(a, d), e = simd_shuffle_down(b, d);
      if (lane + d < width) { a = pair_add(a, c); if (dual) b = pair_add(b, e); }
    }
    if (lane == 0) {
      out[offset] = a.x; out[offset + 1] = a.y;
      if (dual) { out[offset + 2] = b.x; out[offset + 3] = b.y; }
    }
  }
}
kernel void runtime_fill(device uint* dst [[buffer(0)]],
                         constant uint& value [[buffer(1)]],
                         constant uint& n [[buffer(2)]],
                         uint i [[thread_position_in_grid]]) {
  if (i < n) dst[i] = value;
}
)MSL";

class MetalBuffer final : public BufferHandle {
 public:
  MetalBuffer(id<MTLDevice> device, std::size_t bytes) : bytes_(bytes) {
    buffer = [device newBufferWithLength:std::max(bytes, std::size_t{4})
                                options:MTLResourceStorageModeShared];
    KUTTA_REQUIRE(buffer != nil, "Metal buffer allocation failed");
  }
  std::size_t bytes() const override { return bytes_; }
  id<MTLBuffer> buffer;
 private:
  std::size_t bytes_;
};

id<MTLBuffer> metal(const BufferHandle& h) { return static_cast<const MetalBuffer&>(h).buffer; }

class MetalDevice final : public Device {
 public:
  MetalDevice() {
    @autoreleasepool {
      device_ = MTLCreateSystemDefaultDevice();
      KUTTA_REQUIRE(device_ != nil, "no Metal device available");
      queue_ = [device_ newCommandQueue];
      KUTTA_REQUIRE(queue_ != nil, "Metal command queue creation failed");
      MTLCompileOptions* options = [MTLCompileOptions new];
      // Solver kernels need IEEE-accurate sqrt/exp/log/tanh/pow on every OS
      // version; on macOS 15+ mathMode alone still selects fast functions.
      if (@available(macOS 15.0, *)) {
        options.mathMode = MTLMathModeSafe;
        options.mathFloatingPointFunctions = MTLMathFloatingPointFunctionsPrecise;
      } else {
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
        options.fastMathEnabled = NO;
#pragma clang diagnostic pop
      }
      NSError* error = nil;
      std::string joined_source;
      for (std::size_t i = 0; i < metal_kernel_source_part_count; ++i)
        joined_source += metal_kernel_source_parts[i];
      joined_source += reduction_source;
      NSString* source = [NSString stringWithUTF8String:joined_source.c_str()];
      id<MTLLibrary> library = [device_ newLibraryWithSource:source options:options error:&error];
      if (!library) throw Error("Metal compilation: " + std::string(error.localizedDescription.UTF8String));
      auto pipeline = [&](const char* entry) {
        id<MTLFunction> function = [library newFunctionWithName:[NSString stringWithUTF8String:entry]];
        KUTTA_REQUIRE(function != nil, std::string("missing Metal entry: ") + entry);
        id<MTLComputePipelineState> state = [device_ newComputePipelineStateWithFunction:function error:&error];
        if (!state) throw Error("Metal pipeline: " + std::string(error.localizedDescription.UTF8String));
        return state;
      };
      for (std::size_t i = 0; i < kKernelCount; ++i) pipelines_.push_back(pipeline(metal_kernel_entry[i]));
      reduce_ = pipeline("runtime_reduce");
      reduce_into_ = pipeline("runtime_reduce_into");
      fill_ = pipeline("runtime_fill");
    }
  }
  ~MetalDevice() override {
    try { synchronize(); } catch (...) { /* Destruction cannot report GPU failures. */ }
  }
  Backend backend() const override { return Backend::Metal; }
  std::string name() const override { return device_.name.UTF8String; }

  void launch(KernelId id, std::size_t n, std::span<const KernelArg> args) override {
    const auto k = static_cast<std::size_t>(id);
    KUTTA_REQUIRE(k < pipelines_.size(), "launch: bad kernel id");
    check_size(n);
    if (!n) return;
    auto encoder = compute();
    [encoder setComputePipelineState:pipelines_[k]];
    for (std::size_t i = 0; i < args.size(); ++i) {
      const auto& a = args[i];
      if (a.kind == KernelArg::Kind::Buffer) [encoder setBuffer:metal(*a.buffer) offset:0 atIndex:i];
      else [encoder setBytes:a.scalar length:a.scalar_size atIndex:i];
    }
    const auto count = static_cast<std::uint32_t>(n);
    [encoder setBytes:&count length:sizeof count atIndex:args.size()];
    dispatch(pipelines_[k], n);
  }
  void synchronize() override {
    @autoreleasepool {
      commit();
      // Wait for every buffer before reporting, so no GPU work is still writing
      // shared memory when the caller handles the error.
      std::string first_error;
      for (id<MTLCommandBuffer> buffer : pending_) {
        const std::string error = wait(buffer);
        if (first_error.empty()) first_error = error;
      }
      pending_.clear();
      if (!first_error.empty()) throw Error("Metal execution: " + first_error);
    }
  }
  double sum(const Buffer<float>& x, std::size_t n) override {
    KUTTA_REQUIRE(n <= x.size(), "sum: n out of range");
    return reduce(x, x, n, 0);
  }
  double dot(const Buffer<float>& x, const Buffer<float>& y, std::size_t n) override {
    KUTTA_REQUIRE(n <= x.size() && n <= y.size(), "dot: n out of range");
    return reduce(x, y, n, 1);
  }
  double max_abs(const Buffer<float>& x, std::size_t n) override {
    KUTTA_REQUIRE(n <= x.size(), "max_abs: n out of range");
    return reduce(x, x, n, 2);
  }

 protected:
    void reduce_into(const Buffer<float>& x, const Buffer<float>& y,
      const Buffer<float>& u, const Buffer<float>& v, std::size_t n,
      Buffer<float>& dst, std::size_t offset, Reduction mode) override {
    check_reduction(x, y, u, v, n, dst, offset, mode);
    const auto groups = enqueue_reduce(x, y, u, v, n,
        mode == Reduction::Sum ? 0 : (mode == Reduction::Dot ? 1 : 3));
    auto encoder = compute();
    const auto index = static_cast<std::uint32_t>(offset);
    const std::uint32_t dual = mode == Reduction::Dot2 ? 1 : 0;
    [encoder setComputePipelineState:reduce_into_];
    [encoder setBuffer:partial_->buffer offset:0 atIndex:0];
    [encoder setBuffer:metal(*dst.handle()) offset:0 atIndex:1];
    [encoder setBytes:&groups length:sizeof groups atIndex:2];
    [encoder setBytes:&index length:sizeof index atIndex:3];
    [encoder setBytes:&dual length:sizeof dual atIndex:4];
    dispatch(reduce_into_, threads(reduce_into_));
  }
  std::unique_ptr<BufferHandle> alloc_bytes(std::size_t bytes) override {
    return std::make_unique<MetalBuffer>(device_, bytes);
  }
  void upload_bytes(BufferHandle& dst, const void* src, std::size_t bytes) override {
    if (!bytes) return;
    // Shared host writes must wait for prior GPU readers as well as writers.
    synchronize();
    std::memcpy(metal(dst).contents, src, bytes);
  }
  void download_bytes(const BufferHandle& src, void* dst, std::size_t bytes) override {
    synchronize();
    if (bytes) std::memcpy(dst, metal(src).contents, bytes);
  }
  void copy_bytes(BufferHandle& dst, const BufferHandle& src, std::size_t bytes) override {
    if (!bytes || &dst == &src) return;
    @autoreleasepool {
      end_compute();
      open();
      id<MTLBlitCommandEncoder> blit = [command_ blitCommandEncoder];
      [blit copyFromBuffer:metal(src) sourceOffset:0 toBuffer:metal(dst) destinationOffset:0 size:bytes];
      [blit endEncoding];
      count_dispatch();
    }
  }
  void fill_bytes(BufferHandle& dst, const void* value, std::size_t size, std::size_t n) override {
    KUTTA_REQUIRE(size == 4, "Metal fill requires a 32-bit value");
    check_size(n);
    if (!n) return;
    auto encoder = compute();
    auto count = static_cast<std::uint32_t>(n);
    [encoder setComputePipelineState:fill_];
    [encoder setBuffer:metal(dst) offset:0 atIndex:0];
    [encoder setBytes:value length:size atIndex:1];
    [encoder setBytes:&count length:sizeof count atIndex:2];
    dispatch(fill_, n);
  }

 private:
  static void check_size(std::size_t n) {
    KUTTA_REQUIRE(n <= INT_MAX, "Metal launch exceeds signed kernel index range");
  }
  void open() {
    @autoreleasepool {
      if (!command_) {
        command_ = [queue_ commandBuffer];
        KUTTA_REQUIRE(command_ != nil, "Metal command buffer creation failed");
      }
    }
  }
  id<MTLComputeCommandEncoder> compute() {
    @autoreleasepool {
      open();
      if (!encoder_) encoder_ = [command_ computeCommandEncoderWithDispatchType:MTLDispatchTypeSerial];
      KUTTA_REQUIRE(encoder_ != nil, "Metal compute encoder creation failed");
      return encoder_;
    }
  }
  void end_compute() {
    if (encoder_) { [encoder_ endEncoding]; encoder_ = nil; }
  }
  void commit() {
    end_compute();
    if (command_) {
      [command_ commit];
      pending_.push_back(command_);
      command_ = nil;
      dispatches_ = 0;
    }
  }
  static std::string wait(id<MTLCommandBuffer> buffer) {
    [buffer waitUntilCompleted];
    if (buffer.status != MTLCommandBufferStatusError) return {};
    return buffer.error.localizedDescription.UTF8String;
  }
  void count_dispatch() {
    if (++dispatches_ == 256) {
      commit();
      // Bound command-buffer/resource retention for long host-only launch loops,
      // waiting only on the oldest buffer so the GPU stays busy.
      if (pending_.size() >= 8) {
        const std::string error = wait(pending_.front());
        pending_.erase(pending_.begin());
        if (!error.empty()) {
          synchronize();
          throw Error("Metal execution: " + error);
        }
      }
    }
  }
  static NSUInteger threads(id<MTLComputePipelineState> pipeline) {
    return std::min(pipeline.threadExecutionWidth * 8, pipeline.maxTotalThreadsPerThreadgroup);
  }
  void dispatch(id<MTLComputePipelineState> pipeline, std::size_t n) {
    [encoder_ dispatchThreads:MTLSizeMake(n, 1, 1)
       threadsPerThreadgroup:MTLSizeMake(threads(pipeline), 1, 1)];
    count_dispatch();
  }
  std::uint32_t enqueue_reduce(const Buffer<float>& x, const Buffer<float>& y,
      const Buffer<float>& u, const Buffer<float>& v, std::size_t n, std::uint32_t mode) {
    check_size(n);

    const auto count = static_cast<std::uint32_t>(n);
    const auto group_count = std::min<std::size_t>((n + threads(reduce_) - 1) / threads(reduce_), kMaxGroups);
    // Ordered dispatches consume scratch before its next use, including across command buffers.
    if (!partial_) partial_ = std::make_unique<MetalBuffer>(device_, kMaxGroups * 4 * sizeof(float));
    MetalBuffer& partial = *partial_;
    auto encoder = compute();
    [encoder setComputePipelineState:reduce_];
    [encoder setBuffer:metal(*x.handle()) offset:0 atIndex:0];
    [encoder setBuffer:metal(*y.handle()) offset:0 atIndex:1];
    [encoder setBuffer:partial.buffer offset:0 atIndex:2];
    [encoder setBytes:&count length:sizeof count atIndex:3];
    [encoder setBytes:&mode length:sizeof mode atIndex:4];
    [encoder setBuffer:metal(*u.handle()) offset:0 atIndex:5];
    [encoder setBuffer:metal(*v.handle()) offset:0 atIndex:6];
    if (group_count) [encoder dispatchThreadgroups:MTLSizeMake(group_count, 1, 1)
             threadsPerThreadgroup:MTLSizeMake(threads(reduce_), 1, 1)];
    count_dispatch();
    return static_cast<std::uint32_t>(group_count);
  }
  double reduce(const Buffer<float>& x, const Buffer<float>& y, std::size_t n, std::uint32_t mode) {
    const auto group_count = enqueue_reduce(x, y, x, y, n, mode);
    synchronize();
    const float* p = static_cast<const float*>(partial_->buffer.contents);
    double result = 0, correction = 0;
    for (std::size_t i = 0; i < group_count * 2; ++i) {
      if (mode == 2) result = (p[i] > result || std::isnan(p[i])) ? double(p[i]) : result;
      else {
        double v = double(p[i]) - correction;
        double next = result + v;
        // Inf - Inf in the compensation would turn an infinite sum into NaN.
        correction = std::isfinite(next) ? (next - result) - v : 0.0;
        result = next;
      }
    }
    return result;
  }
  id<MTLDevice> device_;
  id<MTLCommandQueue> queue_;
  id<MTLCommandBuffer> command_;
  id<MTLComputeCommandEncoder> encoder_;
  std::vector<id<MTLCommandBuffer>> pending_;
  std::vector<id<MTLComputePipelineState>> pipelines_;
  id<MTLComputePipelineState> reduce_, reduce_into_, fill_;
  unsigned dispatches_ = 0;
  static constexpr std::size_t kMaxGroups = 1024;
  std::unique_ptr<MetalBuffer> partial_;
};
}  // namespace

std::unique_ptr<Device> make_metal_device() { return std::make_unique<MetalDevice>(); }
std::vector<DeviceInfo> enumerate_metal_devices() {
  @autoreleasepool {
    std::vector<DeviceInfo> out;
    id<MTLDevice> preferred = MTLCreateSystemDefaultDevice();
    int ordinal = 0;
    for (id<MTLDevice> d in MTLCopyAllDevices())
      out.push_back({Backend::Metal, d.name.UTF8String,
        std::uint64_t(d.recommendedMaxWorkingSetSize), ordinal++, d.registryID == preferred.registryID});
    return out;
  }
}
}  // namespace kutta
