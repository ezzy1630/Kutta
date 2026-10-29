# CUDA backend

The backend and generated `kernels_cuda.cu` compile cleanly with CUDA 13.0.2 for
sm_75/86/89/120 (Docker `nvidia/cuda:13.0.2-devel-ubuntu24.04` on macOS) but have
never run on an NVIDIA GPU.
Validate with CUDA 12.x/13.x and a C++20 host toolchain (MSVC 2022 or clang/gcc):

```
cmake -S . -B build-cuda -G Ninja -DKUTTA_ENABLE_METAL=OFF -DKUTTA_ENABLE_CUDA=ON
cmake --build build-cuda
./build-cuda/kutta-tests
./build-cuda/bench_bandwidth
```

The default architectures are sm_75/86/89, supported by CUDA 12.x. For Blackwell,
use a supporting toolkit and `-DCMAKE_CUDA_ARCHITECTURES=120`. Tests enumerate
available backends; confirm CUDA actually appears before interpreting success.
All operations use one non-default stream; uploads wait before returning so
caller-owned host storage can be reused. Reductions use two deterministic FP64
device stages. No offline Metal tools are needed by the separate Metal backend.
