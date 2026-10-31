#import <Foundation/Foundation.h>
#import <Metal/Metal.h>

#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include "kernels.hpp"
#include "kutta/core/memory.hpp"
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
// runtime_reduce_into applied to a single partial: lane 0 adds the partial to
// zero, then adds zero once per level of its two shuffle trees. Replaying those
// exact additions keeps the fused one-dispatch result bitwise identical.
float2 finish_single(float2 p, uint width) {
  float2 a = pair_add(float2(0.0f), p);
  for (uint d = width / 2; d > 0; d /= 2) a = pair_add(a, float2(0.0f));
  for (uint d = width / 2; d > 0; d /= 2) a = pair_add(a, float2(0.0f));
  return a;
}
// One threadgroup's (hi, lo) partial, valid in lane 0 of SIMD group 0. Shared
// by runtime_reduce and runtime_sequence so both run identical arithmetic.
void reduce_group(const device float* x, const device float* y, const device float* u,
                  const device float* v, uint n, uint mode, uint gid, uint stride, uint width,
                  uint lane, uint simd, uint threads, threadgroup float2* scratch,
                  threadgroup float2* second, thread float2& a, thread float2& c) {
  a = float2(0.0f); c = float2(0.0f);
  for (uint i = gid; i < n; i += stride) {
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
  }
}
// Writes the group partial(s), or with `fused` (one threadgroup only) the
// finished pair(s) at out[offset].
void reduce_store(float2 a, float2 c, uint mode, uint fused, uint width, uint group,
                  device float2* partial, device float* out, uint offset) {
  if (fused) {
    // Single threadgroup: finish here instead of in runtime_reduce_into.
    a = finish_single(a, width);
    out[offset] = a.x; out[offset + 1] = a.y;
    if (mode == 3) {
      c = finish_single(c, width);
      out[offset + 2] = c.x; out[offset + 3] = c.y;
    }
  } else {
    partial[group] = a;
    if (mode == 3) partial[1024 + group] = c;
  }
}
kernel void runtime_reduce(const device float* x [[buffer(0)]],
                           const device float* y [[buffer(1)]],
                           device float2* partial [[buffer(2)]],
                           constant uint& n [[buffer(3)]],
                           constant uint& mode [[buffer(4)]],
                           const device float* u [[buffer(5)]],
                           const device float* v [[buffer(6)]],
                           device float* out [[buffer(7)]],
                           constant uint& offset [[buffer(8)]],
                           constant uint& fused [[buffer(9)]],
                           uint gid [[thread_position_in_grid]],
                           uint group [[threadgroup_position_in_grid]],
                           uint width [[threads_per_simdgroup]],
                           uint lane [[thread_index_in_simdgroup]],
                           uint simd [[simdgroup_index_in_threadgroup]],
                           uint threads [[threads_per_threadgroup]],
                           uint groups [[threadgroups_per_grid]]) {
  threadgroup float2 scratch[32];
  threadgroup float2 second[32];
  float2 a, c;
  reduce_group(x, y, u, v, n, mode, gid, threads * groups, width, lane, simd, threads,
               scratch, second, a, c);
  if (simd == 0 && lane == 0) reduce_store(a, c, mode, fused, width, group, partial, out, offset);
}
// Replays recorded small operations in one threadgroup, in order, with a
// device-memory barrier between steps: the same visibility that ordered serial
// dispatches give. Kernel bodies are the generated ksrc functions; reductions
// use reduce_group with the layout of a one-group runtime_reduce dispatch.
struct SequenceStep { uint op; uint n; uint mode; uint offset; ulong arg[32]; };
constant uint kSequenceFill = 0xFFFF0001u, kSequenceCopy = 0xFFFF0002u, kSequenceReduce = 0xFFFF0003u;
kernel void runtime_sequence(constant SequenceStep* steps [[buffer(0)]],
                             constant uint& count [[buffer(1)]],
                             uint tid [[thread_index_in_threadgroup]],
                             uint width [[threads_per_simdgroup]],
                             uint lane [[thread_index_in_simdgroup]],
                             uint simd [[simdgroup_index_in_threadgroup]],
                             uint threads [[threads_per_threadgroup]]) {
  threadgroup float2 scratch[32];
  threadgroup float2 second[32];
  for (uint s = 0; s < count; ++s) {
    constant SequenceStep& step = steps[s];
    if (step.op == kSequenceReduce) {
      float2 a, c;
      reduce_group(reinterpret_cast<const device float*>(step.arg[0]),
                   reinterpret_cast<const device float*>(step.arg[1]),
                   reinterpret_cast<const device float*>(step.arg[2]),
                   reinterpret_cast<const device float*>(step.arg[3]),
                   step.n, step.mode, tid, threads, width, lane, simd, threads, scratch, second, a, c);
      if (simd == 0 && lane == 0)
        reduce_store(a, c, step.mode, uint(step.arg[5]), width, 0,
                     reinterpret_cast<device float2*>(step.arg[4]),
                     reinterpret_cast<device float*>(step.arg[4]), step.offset);
    } else if (step.op == kSequenceFill) {
      device uint* dst = reinterpret_cast<device uint*>(step.arg[0]);
      for (uint i = tid; i < step.n; i += threads) dst[i] = uint(step.arg[1]);
    } else if (step.op == kSequenceCopy) {
      device uint* dst = reinterpret_cast<device uint*>(step.arg[0]);
      const device uint* src = reinterpret_cast<const device uint*>(step.arg[1]);
      for (uint i = tid; i < step.n; i += threads) dst[i] = src[i];
    } else {
      for (uint i = tid; i < step.n; i += threads) kutta_sequence_call(step.op, int(i), step.arg);
    }
    threadgroup_barrier(mem_flags::mem_device | mem_flags::mem_threadgroup);
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
    address = buffer.gpuAddress;
  }
  std::size_t bytes() const override { return bytes_; }
  id<MTLBuffer> buffer;
  std::uint64_t address;  // cached: recording a sequence step reads it per argument
 private:
  std::size_t bytes_;
  memory::DeviceBytes tally_{bytes_};
};

id<MTLBuffer> metal(const BufferHandle& h) { return static_cast<const MetalBuffer&>(h).buffer; }
std::uint64_t address(const BufferHandle& h) { return static_cast<const MetalBuffer&>(h).address; }

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
      sequence_ = pipeline("runtime_sequence");
      // Replayed reductions must use runtime_reduce's one-group layout.
      sequence_threads_ = threads(reduce_);
      const char* toggle = std::getenv("KUTTA_METAL_SEQUENCE");
      sequence_enabled_ = !(toggle && std::string(toggle) == "0") &&
          sequence_.maxTotalThreadsPerThreadgroup >= sequence_threads_ &&
          sequence_.threadExecutionWidth == reduce_.threadExecutionWidth;
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
    if (sequence_eligible(n) && args.size() <= 32) {
      SequenceStep step{static_cast<std::uint32_t>(k), static_cast<std::uint32_t>(n), 0, 0, {}};
      for (std::size_t i = 0; i < args.size(); ++i) {
        const auto& a = args[i];
        if (a.kind == KernelArg::Kind::Buffer) {
          step.arg[i] = address(*a.buffer);
          step_buffers_[i] = metal(*a.buffer);
        } else {
          KUTTA_REQUIRE(a.scalar_size == 4, "Metal sequence scalars are 32-bit");
          std::uint32_t bits;
          std::memcpy(&bits, a.scalar, 4);
          step.arg[i] = bits;
          step_buffers_[i] = nil;
        }
      }
      record(step, args.size());
      return;
    }
    flush_sequence();
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
      flush_sequence();
      commit();
      // Wait for every buffer before reporting, so no GPU work is still writing
      // shared memory when the caller handles the error.
      std::string first_error;
      for (id<MTLCommandBuffer> buffer : pending_) {
        const std::string error = wait(buffer);
        if (first_error.empty()) first_error = error;
      }
      pending_.clear();
      reclaim_step_buffers();
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
    const std::uint32_t code = mode == Reduction::Sum ? 0 : (mode == Reduction::Dot ? 1 : 3);
    // Small vectors fit one threadgroup, which finishes the pair itself. This
    // halves the dispatches of the device-scalar Krylov loop on small meshes.
    if (group_count(n) == 1) {
      enqueue_reduce(x, y, u, v, n, code, *partial(), 0, &dst, offset);
      return;
    }
    const auto groups = enqueue_reduce(x, y, u, v, n, code, *partial(), 0);
    auto encoder = compute();
    const auto index = static_cast<std::uint32_t>(offset);
    const std::uint32_t dual = mode == Reduction::Dot2 ? 1 : 0;
    [encoder setComputePipelineState:reduce_into_];
    [encoder setBuffer:partial()->buffer offset:0 atIndex:0];
    [encoder setBuffer:metal(*dst.handle()) offset:0 atIndex:1];
    [encoder setBytes:&groups length:sizeof groups atIndex:2];
    [encoder setBytes:&index length:sizeof index atIndex:3];
    [encoder setBytes:&dual length:sizeof dual atIndex:4];
    dispatch(reduce_into_, threads(reduce_into_));
  }
  void host_reductions(std::span<const HostReduction> requests, std::span<double> out) override {
    // Each request owns a slot of host-finished group partials, so one wait
    // serves the whole batch and every value matches its single-call result.
    constexpr std::size_t slot = kMaxGroups * 2 * sizeof(float);
    // Every batch ends in a wait, so no queued work still references batch_.
    if (!batch_ || batch_->bytes() < requests.size() * slot)
      batch_ = std::make_unique<MetalBuffer>(device_, std::max<std::size_t>(requests.size(), 8) * slot);
    std::vector<std::uint32_t> groups(requests.size()), modes(requests.size());
    for (std::size_t i = 0; i < requests.size(); ++i) {
      const auto& r = requests[i];
      KUTTA_REQUIRE(r.n <= r.x->size() && r.n <= r.y->size(), "reduce_many: n out of range");
      modes[i] = r.kind == HostReduction::Kind::Sum ? 0 : (r.kind == HostReduction::Kind::Dot ? 1 : 2);
      groups[i] = enqueue_reduce(*r.x, *r.y, *r.x, *r.y, r.n, modes[i], *batch_, i * slot);
    }
    synchronize();
    const auto* base = static_cast<const unsigned char*>(batch_->buffer.contents);
    for (std::size_t i = 0; i < requests.size(); ++i)
      out[i] = finish(reinterpret_cast<const float*>(base + i * slot), groups[i], modes[i]);
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
    if (bytes % 4 == 0 && sequence_eligible(bytes / 4)) {
      SequenceStep step{kSequenceCopy, static_cast<std::uint32_t>(bytes / 4), 0, 0, {}};
      step.arg[0] = address(dst); step_buffers_[0] = metal(dst);
      step.arg[1] = address(src); step_buffers_[1] = metal(src);
      record(step, 2);
      return;
    }
    flush_sequence();
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
    if (sequence_eligible(n)) {
      SequenceStep step{kSequenceFill, static_cast<std::uint32_t>(n), 0, 0, {}};
      std::uint32_t bits;
      std::memcpy(&bits, value, 4);
      step.arg[0] = address(dst); step_buffers_[0] = metal(dst);
      step.arg[1] = bits; step_buffers_[1] = nil;
      record(step, 2);
      return;
    }
    flush_sequence();
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
  void count_dispatch(std::size_t weight = 1) {
    dispatches_ += weight;
    if (dispatches_ >= 256) {
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
  // ---- Small-operation sequences ----
  // Tiny meshes are dominated by per-dispatch GPU cost, not by bandwidth.
  // Launches, fills, copies and one-group reductions of at most kSequenceRows
  // elements are recorded and replayed by runtime_sequence in one dispatch.
  // Each step runs the same kernel body or reduction tree as its own dispatch,
  // so results are bitwise identical. Any other operation or wait flushes
  // first, which keeps the Device ordering contract. KUTTA_METAL_SEQUENCE=0
  // disables replay (diagnostics and A/B timing).
  struct SequenceStep { std::uint32_t op, n, mode, offset; std::uint64_t arg[32]; };
  static_assert(sizeof(SequenceStep) == 272, "must match runtime_sequence's SequenceStep");
  static constexpr std::uint32_t kSequenceFill = 0xFFFF0001u, kSequenceCopy = 0xFFFF0002u,
      kSequenceReduce = 0xFFFF0003u;
  static constexpr std::size_t kSequenceRows = 2048, kSequenceSteps = 4096;
  bool sequence_eligible(std::size_t n) const { return sequence_enabled_ && n <= kSequenceRows; }
  // step_buffers_[0, count) holds the step's buffers in argument order (nil
  // for scalars), for residency and for re-encoding a lone step natively.
  void record(const SequenceStep& step, std::size_t count) {
    if (steps_.size() == kSequenceSteps) flush_sequence();
    steps_.push_back(step);
    last_count_ = count;
    // A flush references few distinct buffers; a linear scan beats hashing.
    // resources_ holds the strong references that keep recorded addresses valid.
    for (std::size_t i = 0; i < count; ++i) {
      const void* key = (__bridge const void*)step_buffers_[i];
      if (key && std::find(resource_keys_.begin(), resource_keys_.end(), key) == resource_keys_.end()) {
        resource_keys_.push_back(key);
        resources_.push_back(step_buffers_[i]);
      }
    }
  }
  void flush_sequence() {
    if (steps_.empty()) return;
    @autoreleasepool {
      if (steps_.size() == 1) {
        // A lone step gains nothing from replay; encode its own dispatch.
        // resources_ still holds the only references to buffers whose owners
        // are gone, so clear it only once the encoder has retained them.
        encode_native(steps_.front());
        clear_sequence();
        return;
      }
      auto encoder = compute();
      const auto bytes = steps_.size() * sizeof(SequenceStep);
      id<MTLBuffer> buffer = step_buffer(bytes);
      std::memcpy(buffer.contents, steps_.data(), bytes);
      const auto count = static_cast<std::uint32_t>(steps_.size());
      [encoder setComputePipelineState:sequence_];
      [encoder setBuffer:buffer offset:0 atIndex:0];
      [encoder setBytes:&count length:sizeof count atIndex:1];
      [encoder useResources:resources_.data() count:resources_.size()
                      usage:MTLResourceUsageRead | MTLResourceUsageWrite];
      [encoder dispatchThreadgroups:MTLSizeMake(1, 1, 1)
              threadsPerThreadgroup:MTLSizeMake(sequence_threads_, 1, 1)];
      clear_sequence();
      // Weight by steps so a long replay is committed as early as the
      // dispatches it replaces would have been.
      count_dispatch(count);
    }
  }
  void clear_sequence() {
    steps_.clear();
    resources_.clear();
    resource_keys_.clear();
  }
  void encode_native(const SequenceStep& step) {
    if (step.op == kSequenceCopy) {
      end_compute();
      open();
      id<MTLBlitCommandEncoder> blit = [command_ blitCommandEncoder];
      [blit copyFromBuffer:step_buffers_[1] sourceOffset:0 toBuffer:step_buffers_[0]
         destinationOffset:0 size:std::size_t(step.n) * 4];
      [blit endEncoding];
      count_dispatch();
      return;
    }
    auto encoder = compute();
    if (step.op == kSequenceFill) {
      const auto value = static_cast<std::uint32_t>(step.arg[1]);
      [encoder setComputePipelineState:fill_];
      [encoder setBuffer:step_buffers_[0] offset:0 atIndex:0];
      [encoder setBytes:&value length:sizeof value atIndex:1];
      [encoder setBytes:&step.n length:sizeof step.n atIndex:2];
      dispatch(fill_, step.n);
      return;
    }
    if (step.op == kSequenceReduce) {
      const auto fused = static_cast<std::uint32_t>(step.arg[5]);
      [encoder setComputePipelineState:reduce_];
      for (int i = 0; i < 4; ++i)
        [encoder setBuffer:step_buffers_[i] offset:0 atIndex:i < 2 ? i : i + 3];
      [encoder setBuffer:step_buffers_[4] offset:step.arg[4] - step_buffers_[4].gpuAddress atIndex:2];
      [encoder setBytes:&step.n length:sizeof step.n atIndex:3];
      [encoder setBytes:&step.mode length:sizeof step.mode atIndex:4];
      [encoder setBuffer:step_buffers_[4] offset:0 atIndex:7];
      [encoder setBytes:&step.offset length:sizeof step.offset atIndex:8];
      [encoder setBytes:&fused length:sizeof fused atIndex:9];
      [encoder dispatchThreadgroups:MTLSizeMake(1, 1, 1)
              threadsPerThreadgroup:MTLSizeMake(threads(reduce_), 1, 1)];
      count_dispatch();
      return;
    }
    [encoder setComputePipelineState:pipelines_[step.op]];
    for (std::size_t i = 0; i < last_count_; ++i) {
      if (step_buffers_[i]) [encoder setBuffer:step_buffers_[i] offset:0 atIndex:i];
      else {
        const auto bits = static_cast<std::uint32_t>(step.arg[i]);
        [encoder setBytes:&bits length:sizeof bits atIndex:i];
      }
    }
    [encoder setBytes:&step.n length:sizeof step.n atIndex:last_count_];
    dispatch(pipelines_[step.op], step.n);
  }
  // Step lists live in shared buffers until their command buffer completes.
  id<MTLBuffer> step_buffer(std::size_t bytes) {
    reclaim_step_buffers();
    constexpr std::size_t size = kSequenceSteps * sizeof(SequenceStep);
    KUTTA_REQUIRE(bytes <= size, "Metal sequence exceeds its step buffer");
    id<MTLBuffer> buffer = nil;
    if (!free_step_buffers_.empty()) {
      buffer = free_step_buffers_.back();
      free_step_buffers_.pop_back();
    } else {
      buffer = [device_ newBufferWithLength:size options:MTLResourceStorageModeShared];
      KUTTA_REQUIRE(buffer != nil, "Metal sequence buffer allocation failed");
    }
    used_step_buffers_.push_back({buffer, command_});
    return buffer;
  }
  void reclaim_step_buffers() {
    std::erase_if(used_step_buffers_, [&](const auto& used) {
      if (used.second.status < MTLCommandBufferStatusCompleted) return false;
      free_step_buffers_.push_back(used.first);
      return true;
    });
  }
  std::size_t group_count(std::size_t n) {
    return std::min<std::size_t>((n + threads(reduce_) - 1) / threads(reduce_), kMaxGroups);
  }
  // Ordered dispatches consume scratch before its next use, including across command buffers.
  MetalBuffer* partial() {
    if (!partial_) partial_ = std::make_unique<MetalBuffer>(device_, kMaxGroups * 4 * sizeof(float));
    return partial_.get();
  }
  // Writes per-group pairs to `partial` at byte `partial_offset`, or with
  // `out` set (one group only) the finished pair(s) to out[offset].
  std::uint32_t enqueue_reduce(const Buffer<float>& x, const Buffer<float>& y,
      const Buffer<float>& u, const Buffer<float>& v, std::size_t n, std::uint32_t mode,
      MetalBuffer& partial, std::size_t partial_offset,
      const Buffer<float>* out = nullptr, std::size_t offset = 0) {
    check_size(n);
    const auto count = static_cast<std::uint32_t>(n);
    const auto groups = group_count(n);
    KUTTA_REQUIRE(!out || groups == 1, "Metal fused reduction needs one threadgroup");
    const auto out_offset = static_cast<std::uint32_t>(offset);
    const std::uint32_t fused = out ? 1 : 0;
    if (groups == 1 && sequence_eligible(n)) {
      SequenceStep step{kSequenceReduce, count, mode, out_offset, {}};
      const Buffer<float>* inputs[]{&x, &y, &u, &v};
      for (int i = 0; i < 4; ++i) {
        step.arg[i] = address(*inputs[i]->handle());
        step_buffers_[i] = metal(*inputs[i]->handle());
      }
      step_buffers_[4] = out ? metal(*out->handle()) : partial.buffer;
      step.arg[4] = out ? address(*out->handle()) : partial.address + partial_offset;
      step.arg[5] = fused; step_buffers_[5] = nil;
      record(step, 6);
      return 1;
    }
    flush_sequence();
    auto encoder = compute();
    [encoder setComputePipelineState:reduce_];
    [encoder setBuffer:metal(*x.handle()) offset:0 atIndex:0];
    [encoder setBuffer:metal(*y.handle()) offset:0 atIndex:1];
    [encoder setBuffer:partial.buffer offset:partial_offset atIndex:2];
    [encoder setBytes:&count length:sizeof count atIndex:3];
    [encoder setBytes:&mode length:sizeof mode atIndex:4];
    [encoder setBuffer:metal(*u.handle()) offset:0 atIndex:5];
    [encoder setBuffer:metal(*v.handle()) offset:0 atIndex:6];
    [encoder setBuffer:out ? metal(*out->handle()) : partial.buffer offset:0 atIndex:7];
    [encoder setBytes:&out_offset length:sizeof out_offset atIndex:8];
    [encoder setBytes:&fused length:sizeof fused atIndex:9];
    if (groups) [encoder dispatchThreadgroups:MTLSizeMake(groups, 1, 1)
             threadsPerThreadgroup:MTLSizeMake(threads(reduce_), 1, 1)];
    count_dispatch();
    return static_cast<std::uint32_t>(groups);
  }
  static double finish(const float* p, std::size_t group_count, std::uint32_t mode) {
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
  double reduce(const Buffer<float>& x, const Buffer<float>& y, std::size_t n, std::uint32_t mode) {
    const auto groups = enqueue_reduce(x, y, x, y, n, mode, *partial(), 0);
    synchronize();
    return finish(static_cast<const float*>(partial_->buffer.contents), groups, mode);
  }
  id<MTLDevice> device_;
  id<MTLCommandQueue> queue_;
  id<MTLCommandBuffer> command_;
  id<MTLComputeCommandEncoder> encoder_;
  std::vector<id<MTLCommandBuffer>> pending_;
  std::vector<id<MTLComputePipelineState>> pipelines_;
  id<MTLComputePipelineState> reduce_, reduce_into_, fill_, sequence_;
  std::size_t dispatches_ = 0;
  bool sequence_enabled_ = false;
  NSUInteger sequence_threads_ = 0;
  std::vector<SequenceStep> steps_;
  std::size_t last_count_ = 0;
  // Unretained: every recorded buffer is also strongly held in resources_.
  __unsafe_unretained id<MTLBuffer> step_buffers_[32];
  std::vector<id<MTLBuffer>> resources_;
  std::vector<const void*> resource_keys_;
  std::vector<std::pair<id<MTLBuffer>, id<MTLCommandBuffer>>> used_step_buffers_;
  std::vector<id<MTLBuffer>> free_step_buffers_;
  static constexpr std::size_t kMaxGroups = 1024;
  std::unique_ptr<MetalBuffer> partial_, batch_;
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
