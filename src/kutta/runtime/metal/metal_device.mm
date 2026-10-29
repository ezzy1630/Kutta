#import <Foundation/Foundation.h>
#import <Metal/Metal.h>

#include <algorithm>
#include <climits>
#include <cstring>
#include <vector>

#include "kernels.hpp"
#include "kutta/runtime/metal/metal_kernel_table.hpp"

namespace kutta {
namespace {

// TwoSum keeps the rounding residual at every addition, including across SIMD
// lanes and threadgroups. Dot products also retain the FMA product residual.
// With fast math disabled these float pairs have roughly 48 bits of precision;
// the host finishes the bounded set of (hi, lo) partials in double.
constexpr const char* reduction_source = R"MSL(
#include <metal_stdlib>
using namespace metal;
float2 pair_add(float2 a, float2 b) {
  float s = a.x + b.x;
  float v = s - a.x;
  float e = (a.x - (s - v)) + (b.x - v);
  e = e + a.y + b.y;
  float h = s + e;
  return float2(h, e - (h - s));
}
kernel void runtime_reduce(const device float* x [[buffer(0)]],
                           const device float* y [[buffer(1)]],
                           device float2* partial [[buffer(2)]],
                           constant uint& n [[buffer(3)]],
                           constant uint& mode [[buffer(4)]],
                           uint gid [[thread_position_in_grid]],
                           uint tid [[thread_index_in_threadgroup]],
                           uint group [[threadgroup_position_in_grid]],
                           uint width [[threads_per_simdgroup]],
                           uint lane [[thread_index_in_simdgroup]],
                           uint simd [[simdgroup_index_in_threadgroup]],
                           uint threads [[threads_per_threadgroup]],
                           uint groups [[threadgroups_per_grid]]) {
  threadgroup float2 scratch[32];
  float2 a = float2(0.0f);
  for (uint i = gid; i < n; i += threads * groups) {
    if (mode == 2) a.x = max(a.x, abs(x[i]));
    else {
      float h = mode == 1 ? x[i] * y[i] : x[i];
      float l = mode == 1 ? fma(x[i], y[i], -h) : 0.0f;
      a = pair_add(a, float2(h, l));
    }
  }
  for (uint d = width / 2; d > 0; d /= 2) {
    float2 b = simd_shuffle_down(a, d);
    if (lane + d < width) a = mode == 2 ? float2(max(a.x, b.x), 0) : pair_add(a, b);
  }
  if (lane == 0) scratch[simd] = a;
  threadgroup_barrier(mem_flags::mem_threadgroup);
  if (simd == 0) {
    a = lane < threads / width ? scratch[lane] : float2(0.0f);
    for (uint d = width / 2; d > 0; d /= 2) {
      float2 b = simd_shuffle_down(a, d);
      if (lane + d < width) a = mode == 2 ? float2(max(a.x, b.x), 0) : pair_add(a, b);
    }
    if (lane == 0) partial[group] = a;
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
      if (@available(macOS 15.0, *)) options.mathMode = MTLMathModeSafe;
      else {
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
        options.fastMathEnabled = NO;
#pragma clang diagnostic pop
      }
      NSError* error = nil;
      NSString* source = [[NSString stringWithUTF8String:metal_kernel_source]
          stringByAppendingString:[NSString stringWithUTF8String:reduction_source]];
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
      for (id<MTLCommandBuffer> buffer : pending_) {
        [buffer waitUntilCompleted];
        if (buffer.status == MTLCommandBufferStatusError) {
          const std::string message = buffer.error.localizedDescription.UTF8String;
          pending_.clear();
          throw Error("Metal execution: " + message);
        }
      }
      pending_.clear();
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
  void count_dispatch() {
    if (++dispatches_ == 256) {
      commit();
      // Bound command-buffer/resource retention for long host-only launch loops.
      if (pending_.size() >= 8) synchronize();
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
  double reduce(const Buffer<float>& x, const Buffer<float>& y, std::size_t n, std::uint32_t mode) {
    check_size(n);
    if (!n) { synchronize(); return 0; }
    const auto count = static_cast<std::uint32_t>(n);
    const auto group_count = std::min<std::size_t>((n + threads(reduce_) - 1) / threads(reduce_), 1024);
    MetalBuffer partial(device_, group_count * 2 * sizeof(float));
    auto encoder = compute();
    [encoder setComputePipelineState:reduce_];
    [encoder setBuffer:metal(*x.handle()) offset:0 atIndex:0];
    [encoder setBuffer:metal(*y.handle()) offset:0 atIndex:1];
    [encoder setBuffer:partial.buffer offset:0 atIndex:2];
    [encoder setBytes:&count length:sizeof count atIndex:3];
    [encoder setBytes:&mode length:sizeof mode atIndex:4];
    [encoder dispatchThreadgroups:MTLSizeMake(group_count, 1, 1)
             threadsPerThreadgroup:MTLSizeMake(threads(reduce_), 1, 1)];
    count_dispatch();
    synchronize();
    const float* p = static_cast<const float*>(partial.buffer.contents);
    double result = 0, correction = 0;
    for (std::size_t i = 0; i < group_count * 2; ++i) {
      if (mode == 2) result = std::max(result, double(p[i]));
      else {
        double v = double(p[i]) - correction;
        double next = result + v;
        correction = (next - result) - v;
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
  id<MTLComputePipelineState> reduce_, fill_;
  unsigned dispatches_ = 0;
};
}  // namespace

std::unique_ptr<Device> make_metal_device() { return std::make_unique<MetalDevice>(); }
}  // namespace kutta
