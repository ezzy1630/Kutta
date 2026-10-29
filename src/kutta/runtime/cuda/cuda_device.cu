#include <algorithm>
#include <climits>
#include <cstring>

#include "kernels.hpp"
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
  ~CudaBuffer() override { cudaFree(data); }
  std::size_t bytes() const override { return bytes_; }
  void* data = nullptr;
 private:
  std::size_t bytes_;
};

__global__ void fill_words(unsigned int* dst, unsigned int value, std::size_t n) {
  const std::size_t i = std::size_t(blockIdx.x) * blockDim.x + threadIdx.x;
  if (i < n) dst[i] = value;
}

// Both stages use a fixed tree and no atomics; scheduling cannot change the sum.
__global__ void reduce_first(const float* x, const float* y, std::size_t n,
                             double* partial, int mode) {
  __shared__ double values[kCudaBlock];
  double a = 0;
  for (std::size_t i = std::size_t(blockIdx.x) * blockDim.x + threadIdx.x;
       i < n; i += std::size_t(blockDim.x) * gridDim.x) {
    if (mode == 2) a = fmax(a, fabs(double(x[i])));
    else a += mode == 1 ? double(x[i]) * double(y[i]) : double(x[i]);
  }
  values[threadIdx.x] = a;
  __syncthreads();
  for (unsigned int d = blockDim.x / 2; d; d /= 2) {
    if (threadIdx.x < d) values[threadIdx.x] = mode == 2
        ? fmax(values[threadIdx.x], values[threadIdx.x + d])
        : values[threadIdx.x] + values[threadIdx.x + d];
    __syncthreads();
  }
  if (!threadIdx.x) partial[blockIdx.x] = values[0];
}

__global__ void reduce_last(const double* partial, unsigned int n, double* result, int mode) {
  __shared__ double values[kCudaBlock];
  double a = 0;
  for (unsigned int i = threadIdx.x; i < n; i += blockDim.x)
    a = mode == 2 ? fmax(a, partial[i]) : a + partial[i];
  values[threadIdx.x] = a;
  __syncthreads();
  for (unsigned int d = blockDim.x / 2; d; d /= 2) {
    if (threadIdx.x < d) values[threadIdx.x] = mode == 2
        ? fmax(values[threadIdx.x], values[threadIdx.x + d])
        : values[threadIdx.x] + values[threadIdx.x + d];
    __syncthreads();
  }
  if (!threadIdx.x) *result = values[0];
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
  double reduce(const Buffer<float>& x, const Buffer<float>& y, std::size_t n, int mode) {
    activate();
    if (!n) { synchronize(); return 0; }
    const auto blocks = static_cast<unsigned int>(std::min<std::size_t>((n + kCudaBlock - 1) / kCudaBlock, 1024));
    CudaBuffer partial((blocks + 1) * sizeof(double));
    auto* p = static_cast<double*>(partial.data);
    reduce_first<<<blocks, kCudaBlock, 0, stream_>>>(
        static_cast<const float*>(cuda_ptr(x.handle())),
        static_cast<const float*>(cuda_ptr(y.handle())), n, p, mode);
    check(cudaGetLastError());
    reduce_last<<<1, kCudaBlock, 0, stream_>>>(p, blocks, p + blocks, mode);
    check(cudaGetLastError());
    double result;
    check(cudaMemcpyAsync(&result, p + blocks, sizeof result, cudaMemcpyDeviceToHost, stream_));
    synchronize();
    return result;
  }
  int ordinal_ = 0;
  cudaStream_t stream_ = nullptr;
  std::string name_;
};
}  // namespace

void* cuda_ptr(BufferHandle* h) { return static_cast<CudaBuffer*>(h)->data; }
std::unique_ptr<Device> make_cuda_device() { return std::make_unique<CudaDevice>(); }
}  // namespace kutta
