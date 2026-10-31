#include <algorithm>
#include <climits>
#include <cstring>

#include "kernels.hpp"
#include "kutta/core/memory.hpp"
#include "kutta/runtime/cuda/cuda_kernel_table.hpp"

namespace kutta {
namespace {

void check(cudaError_t status) {
  if (status != cudaSuccess) throw Error(std::string("CUDA: ") + cudaGetErrorString(status));
}

class CudaBuffer final : public BufferHandle {
 public:
  explicit CudaBuffer(std::size_t bytes) : bytes_(bytes) {
    check(cudaMalloc(&data, std::max(bytes, std::size_t{4})));
  }
  // cudaFree synchronizes the device, so freeing a buffer still referenced by
  // queued work is safe (and is why reduction scratch is reused, not freed).
  ~CudaBuffer() override { cudaFree(data); }
  std::size_t bytes() const override { return bytes_; }
  void* data = nullptr;
 private:
  std::size_t bytes_;
  memory::DeviceBytes tally_{bytes_};
};

__global__ void fill_words(unsigned int* dst, unsigned int value, std::size_t n) {
  const std::size_t i = std::size_t(blockIdx.x) * blockDim.x + threadIdx.x;
  if (i < n) dst[i] = value;
}

// NaN-propagating max (fmax drops NaN), so a diverged field never reports a finite norm.
__device__ double nan_max(double a, double b) { return (b > a || isnan(b)) ? b : a; }

// Both stages use a fixed tree and no atomics; scheduling cannot change the sum.
__global__ void reduce_first(const float* x, const float* y, std::size_t n,
                             double* partial, int mode, const float* u, const float* v) {
  __shared__ double values[kCudaBlock], second[kCudaBlock];
  double a = 0, b = 0;
  for (std::size_t i = std::size_t(blockIdx.x) * blockDim.x + threadIdx.x;
       i < n; i += std::size_t(blockDim.x) * gridDim.x) {
    if (mode == 2) a = nan_max(a, fabs(double(x[i])));
    else a += mode == 1 || mode == 3 ? double(x[i]) * double(y[i]) : double(x[i]);
    if (mode == 3) b += double(u[i]) * double(v[i]);
  }
  values[threadIdx.x] = a;
  second[threadIdx.x] = b;
  __syncthreads();
  for (unsigned int d = blockDim.x / 2; d; d /= 2) {
    if (threadIdx.x < d) values[threadIdx.x] = mode == 2
        ? nan_max(values[threadIdx.x], values[threadIdx.x + d])
        : values[threadIdx.x] + values[threadIdx.x + d];
    if (mode == 3 && threadIdx.x < d) second[threadIdx.x] += second[threadIdx.x + d];
    __syncthreads();
  }
  if (!threadIdx.x) {
    partial[blockIdx.x] = values[0];
    if (mode == 3) partial[1024 + blockIdx.x] = second[0];
  }
}

__global__ void reduce_last(const double* partial, unsigned int n, double* result, int mode) {
  __shared__ double values[kCudaBlock];
  double a = 0;
  for (unsigned int i = threadIdx.x; i < n; i += blockDim.x)
    a = mode == 2 ? nan_max(a, partial[i]) : a + partial[i];
  values[threadIdx.x] = a;
  __syncthreads();
  for (unsigned int d = blockDim.x / 2; d; d /= 2) {
    if (threadIdx.x < d) values[threadIdx.x] = mode == 2
        ? nan_max(values[threadIdx.x], values[threadIdx.x + d])
        : values[threadIdx.x] + values[threadIdx.x + d];
    __syncthreads();
  }
  if (!threadIdx.x) *result = values[0];
}

__global__ void reduce_pair(const double* partial, unsigned int n, float* out,
                            std::size_t offset, bool dual) {
  // Fixed tree, FP64 until the final pair encoding.
  __shared__ double first[kCudaBlock], second[kCudaBlock];
  double a = 0, b = 0;
  for (unsigned int j = threadIdx.x; j < n; j += blockDim.x) {
    a += partial[j];
    if (dual) b += partial[1024 + j];
  }
  first[threadIdx.x] = a; second[threadIdx.x] = b;
  __syncthreads();
  for (unsigned int d = blockDim.x / 2; d; d /= 2) {
    if (threadIdx.x < d) {
      first[threadIdx.x] += first[threadIdx.x + d];
      if (dual) second[threadIdx.x] += second[threadIdx.x + d];
    }
    __syncthreads();
  }
  if (threadIdx.x != 0) return;
  a = first[0]; b = second[0];
  const float hi = static_cast<float>(a);
  out[offset] = hi;
  out[offset + 1] = isfinite(hi) ? static_cast<float>(a - double(hi)) : 0;
  if (dual) {
    const float hi2 = static_cast<float>(b);
    out[offset + 2] = hi2;
    out[offset + 3] = isfinite(hi2) ? static_cast<float>(b - double(hi2)) : 0;
  }
}

class CudaDevice final : public Device {
 public:
  CudaDevice() {
    int count = 0;
    check(cudaGetDeviceCount(&count));
    KUTTA_REQUIRE(count > 0, "no CUDA device available");
    check(cudaGetDevice(&ordinal_));
    cudaDeviceProp properties{};
    check(cudaGetDeviceProperties(&properties, ordinal_));
    name_ = properties.name;
    check(cudaStreamCreateWithFlags(&stream_, cudaStreamNonBlocking));
  }
  ~CudaDevice() override {
    cudaSetDevice(ordinal_);
    cudaStreamSynchronize(stream_);
    cudaStreamDestroy(stream_);
  }
  Backend backend() const override { return Backend::Cuda; }
  std::string name() const override { return name_; }
  void synchronize() override { activate(); check(cudaStreamSynchronize(stream_)); }
  void launch(KernelId id, std::size_t n, std::span<const KernelArg> args) override {
    const auto k = static_cast<std::size_t>(id);
    KUTTA_REQUIRE(k < kKernelCount, "launch: bad kernel id");
    KUTTA_REQUIRE(n <= INT_MAX, "CUDA launch exceeds signed kernel index range");
    if (!n) return;
    activate();
    cuda_kernel_table[k](stream_, n, args.data());
    check(cudaGetLastError());
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
    const auto blocks = enqueue_reduce(x, y, u, v, n,
        mode == Reduction::Sum ? 0 : (mode == Reduction::Dot ? 1 : 3));
    reduce_pair<<<1, kCudaBlock, 0, stream_>>>(static_cast<double*>(partial_->data), blocks,
        static_cast<float*>(cuda_ptr(dst.handle())), offset, mode == Reduction::Dot2);
    check(cudaGetLastError());
  }
  void host_reductions(std::span<const HostReduction> requests, std::span<double> out) override {
    // Each finished value gets its own result slot; one copy and one stream
    // wait serve up to kMaxBatch values, each identical to its single call.
    for (std::size_t first = 0; first < requests.size(); first += kMaxBatch) {
      const auto count = std::min(kMaxBatch, requests.size() - first);
      for (std::size_t i = 0; i < count; ++i) {
        const auto& r = requests[first + i];
        KUTTA_REQUIRE(r.n <= r.x->size() && r.n <= r.y->size(), "reduce_many: n out of range");
        const int mode = r.kind == HostReduction::Kind::Sum ? 0 : (r.kind == HostReduction::Kind::Dot ? 1 : 2);
        const auto blocks = enqueue_reduce(*r.x, *r.y, *r.x, *r.y, r.n, mode);
        auto* p = static_cast<double*>(partial_->data);
        reduce_last<<<1, kCudaBlock, 0, stream_>>>(p, blocks, p + 2 * kMaxBlocks + i, mode);
        check(cudaGetLastError());
      }
      check(cudaMemcpyAsync(out.data() + first, static_cast<double*>(partial_->data) + 2 * kMaxBlocks,
                            count * sizeof(double), cudaMemcpyDeviceToHost, stream_));
      synchronize();
    }
  }
  std::unique_ptr<BufferHandle> alloc_bytes(std::size_t bytes) override {
    activate();
    return std::make_unique<CudaBuffer>(bytes);
  }
  void upload_bytes(BufferHandle& dst, const void* src, std::size_t bytes) override {
    activate();
    if (bytes) check(cudaMemcpyAsync(cuda_ptr(&dst), src, bytes, cudaMemcpyHostToDevice, stream_));
    // Caller-owned host memory may be released or overwritten on return.
    synchronize();
  }
  void download_bytes(const BufferHandle& src, void* dst, std::size_t bytes) override {
    activate();
    if (bytes) check(cudaMemcpyAsync(dst, static_cast<const CudaBuffer&>(src).data,
                                    bytes, cudaMemcpyDeviceToHost, stream_));
    synchronize();
  }
  void copy_bytes(BufferHandle& dst, const BufferHandle& src, std::size_t bytes) override {
    activate();
    if (bytes && &dst != &src) check(cudaMemcpyAsync(cuda_ptr(&dst),
        static_cast<const CudaBuffer&>(src).data, bytes, cudaMemcpyDeviceToDevice, stream_));
  }
  void fill_bytes(BufferHandle& dst, const void* value, std::size_t size, std::size_t n) override {
    KUTTA_REQUIRE(size == 4 && n <= INT_MAX, "CUDA fill requires 32-bit values and indices");
    if (!n) return;
    activate();
    unsigned int word;
    std::memcpy(&word, value, sizeof word);
    const auto blocks = static_cast<unsigned int>((n + kCudaBlock - 1) / kCudaBlock);
    fill_words<<<blocks, kCudaBlock, 0, stream_>>>(static_cast<unsigned int*>(cuda_ptr(&dst)), word, n);
    check(cudaGetLastError());
  }
 private:
  void activate() const { check(cudaSetDevice(ordinal_)); }
  unsigned int enqueue_reduce(const Buffer<float>& x, const Buffer<float>& y,
      const Buffer<float>& u, const Buffer<float>& v, std::size_t n, int mode) {
    activate();

    const auto blocks = static_cast<unsigned int>(std::min<std::size_t>((n + kCudaBlock - 1) / kCudaBlock, kMaxBlocks));
    // Scratch reuse is ordered on stream_, including the final pair write.
    if (!partial_) partial_ = std::make_unique<CudaBuffer>((2 * kMaxBlocks + kMaxBatch) * sizeof(double));
    auto* p = static_cast<double*>(partial_->data);
    if (blocks) reduce_first<<<blocks, kCudaBlock, 0, stream_>>>(
        static_cast<const float*>(cuda_ptr(x.handle())),
        static_cast<const float*>(cuda_ptr(y.handle())), n, p, mode,
        static_cast<const float*>(cuda_ptr(u.handle())), static_cast<const float*>(cuda_ptr(v.handle())));
    check(cudaGetLastError());
    return blocks;
  }
  double reduce(const Buffer<float>& x, const Buffer<float>& y, std::size_t n, int mode) {
    const auto blocks = enqueue_reduce(x, y, x, y, n, mode);
    auto* p = static_cast<double*>(partial_->data);
    reduce_last<<<1, kCudaBlock, 0, stream_>>>(p, blocks, p + 2 * kMaxBlocks, mode);
    check(cudaGetLastError());
    double result;
    check(cudaMemcpyAsync(&result, p + 2 * kMaxBlocks, sizeof result, cudaMemcpyDeviceToHost, stream_));
    synchronize();
    return result;
  }
  static constexpr std::size_t kMaxBlocks = 1024, kMaxBatch = 16;
  std::unique_ptr<CudaBuffer> partial_;
  int ordinal_ = 0;
  cudaStream_t stream_ = nullptr;
  std::string name_;
};
}  // namespace

void* cuda_ptr(BufferHandle* h) { return static_cast<CudaBuffer*>(h)->data; }
std::unique_ptr<Device> make_cuda_device() { return std::make_unique<CudaDevice>(); }
std::vector<DeviceInfo> enumerate_cuda_devices() {
  std::vector<DeviceInfo> out;
  int count = 0, preferred = 0;
  if (cudaGetDeviceCount(&count) != cudaSuccess) return out;
  if (cudaGetDevice(&preferred) != cudaSuccess) preferred = 0;
  for (int i = 0; i < count; ++i) {
    cudaDeviceProp p{};
    if (cudaGetDeviceProperties(&p, i) == cudaSuccess)
      out.push_back({Backend::Cuda, p.name, std::uint64_t(p.totalGlobalMem), i, i == preferred});
  }
  return out;
}
}  // namespace kutta
