# NVIDIA build, test and benchmark path for Kutta

Research date: 2026-10-29. Prices change weekly on the marketplaces; every figure below has a source link and a confidence label. "Verified" means read from the vendor's own page or tested on this Mac. "Tracker" means a third-party price aggregator. "Unverified" means the number came from secondary sources, or is my estimate.

## Bottom line

- **(a) Compile checks:** Use two checks, both free.
  - **Local, on every change:** run nvcc in an arm64 Linux Docker container on the Mac. Device code compiles the same for every `sm_*` target whatever the host CPU.
  - **On every push:** run GitHub Actions with `ubuntu-24.04` and `windows-2022` jobs. The Windows job uses MSVC 2022 with CUDA installed by `Jimver/cuda-toolkit`. Compiling needs no GPU.
- **(b) Correctness on a real NVIDIA GPU:** Use Modal Functions on an **L4** (Ada, sm_89, the same architecture as RTX 40). The agent drives it from the Mac with `modal run`. About 20 L4-hours a month fit inside Modal's **$30/month free credit**.
- **(c) Benchmarks:**
  - **Ada (RTX 40):** Modal L4 (300 GB/s, close to a 4060/4070 laptop) and L40S (864 GB/s, close to a 4090/5080/5090 laptop).
  - **Blackwell sm_120 (RTX 50):** rent an **RTX 5070 Ti on Vast.ai**, about $0.13–0.24/hr. It uses the same GB203 die and the same 896 GB/s as the RTX 5080 laptop GPU. The fallback is an RTX 5090 on RunPod at $0.69/hr.
- **Toolkit:** pin **CUDA 13.0.x** (13.0.2 tested here). Ship SASS for `sm_75, sm_86, sm_89, sm_120` plus `compute_120` PTX. The user's machine needs NVIDIA driver **R580 or newer**.
- **Expected spend:** about **$0–5/month** for ~20 GPU-hours. That is ~15 h on Modal L4, inside the free credit, plus ~5 h of RTX 5070 Ti on Vast.ai at about $1. The only cash outlay is a prepaid Vast.ai top-up, which needs the owner's approval.

## 1. Target hardware: what "representative" means

Kutta's solver kernels are stencil/flux sweeps, so they are expected to be memory-bandwidth-bound. Peak DRAM bandwidth is therefore the first-order axis for comparing GPUs. Laptop GPUs also run at 60–175 W, so sustained clocks vary by chassis. Bandwidth-bound kernels are hurt less by this than compute-bound ones.

### User-class GPUs (laptop unless noted)

| GPU | Arch / SM | VRAM | Peak BW (GB/s) | Source |
|---|---|---|---|---|
| RTX 4060 / 4070 Laptop | Ada AD107/AD106, sm_89 | 8 GB | 256 (128-bit × 16 Gbps) | [Notebookcheck 4070L](https://www.notebookcheck.net/NVIDIA-GeForce-RTX-4090-Laptop-GPU-vs-NVIDIA-GeForce-RTX-4070-Laptop-GPU_11437_11453.247598.0.html) (bus width + clock; BW is computed) |
| RTX 4080 Laptop | Ada AD104, sm_89 | 12 GB | 432 | [gpu-monkey](https://www.gpu-monkey.com/en/compare_gpu-nvidia_geforce_rtx_4080_mobile_laptop_gpu-vs-nvidia_geforce_rtx_4090_mobile_laptop_gpu) |
| RTX 4090 Laptop | Ada AD103, sm_89 | 16 GB | 576 (NVIDIA spec; Notebookcheck lists 20 Gbps ≈ 640) | [Notebookcheck](https://www.notebookcheck.net/NVIDIA-GeForce-RTX-4090-Laptop-GPU-Benchmarks-and-Specs.675091.0.html) |
| RTX 5080 Laptop | Blackwell GB203, sm_120 | 16 GB | 896 | [Notebookcheck](https://www.notebookcheck.net/GeForce-RTX-5090-Laptop-vs-GeForce-RTX-5080-Laptop_12928_12927.247598.0.html) |
| RTX 5090 Laptop | Blackwell GB203, sm_120 | 24 GB | 896 | same |
| RTX 4090 (desktop) | Ada AD102, sm_89 | 24 GB | 1008 | NVIDIA spec |
| RTX 5090 (desktop) | Blackwell GB202, sm_120 | 32 GB | 1792 | [Spheron](https://www.spheron.network/blog/nvidia-rtx-5090-specs/) |

### Rentable GPUs and how close they are to the user's

| Rentable GPU | Arch / SM | Peak BW (GB/s) | Closest user GPU |
|---|---|---|---|
| T4 | Turing, sm_75 | 320 | Oldest supported arch (RTX 20). Useful for coverage, not for perf |
| L4 | Ada AD104, sm_89 | 300 | **RTX 4060/4070 Laptop** (same arch, ~1.2× BW) |
| A10 / A10G | Ampere GA102, sm_86 | 600 | RTX 30 arch; BW like a 4090 Laptop, but older arch |
| L40S | Ada AD102, sm_89 | 864 | **RTX 4090 Laptop to 5080/5090 Laptop BW class**, Ada arch |
| RTX 4090 | Ada, sm_89 | 1008 | Upper bound for Ada; ~1.75× a 4090 Laptop |
| RTX 5060 Ti (desktop) | Blackwell GB206, sm_120 | 448 | Low-end RTX 50 |
| **RTX 5070 Ti (desktop)** | Blackwell **GB203**, sm_120 | **896** | **RTX 5080/5090 Laptop: same die, same BW** |
| RTX 5090 | Blackwell GB202, sm_120 | 1792 | 2× a 5090 Laptop |
| RTX PRO 6000 (Modal) | Blackwell GB202, sm_120 | ~1600 (Server Ed.) – 1792 | 2× a 5090 Laptop; expensive |

T4, L4, A10 and L40S bandwidths are NVIDIA datasheet figures (unverified on this pass; they are long-stable published values). The RTX PRO 6000 Server Edition figure is from [Flopper](https://flopper.io/compare/nvidia-rtx-pro-6000-blackwell-server-edition-vs-nvidia-geforce-rtx-5090-32gb).

### Reporting bandwidth honestly

Every benchmark run should do three things:

1. Measure achieved DRAM bandwidth on the same box with a STREAM-style copy/triad kernel built into the benchmark binary.
2. Report solver throughput as cells/s and as a fraction of that *measured* bandwidth.
3. Label anything about the user's laptop as a projection: `measured_fraction × laptop_peak_BW`. Never present a desktop or data-centre number as a laptop number.

The RTX 5070 Ti and L4 cases need only small extrapolation; an RTX 5090 or L40S run needs explicit scaling.

## 2. Provider comparison

Billing granularity, Windows availability and driving method matter as much as $/hr. All prices are USD/hour for one GPU.

### Modal

- **GPUs:** T4 $0.59, L4 $0.80, A10 $1.10, L40S $1.95, RTX PRO 6000 $3.03, plus A100/H100/B200. Per-second billing ([pricing, verified](https://modal.com/pricing)).
- **Other compute:** CPU $0.047/core-h and memory $0.008/GiB-h are billed on top.
- **Free tier:** Starter plan, $0 plus compute, includes **$30/month free credits**.
- **Multipliers:** region pinning costs 1.15–1.75×. "Non-preemptible" costs 3×, and Sandboxes are billed at that rate for CPU/memory (GPU stays at the standard rate).
- **Egress:** $0.04/GiB.
- **Windows:** no. Linux containers only.
- **Driver:** host driver **610.57.04, CUDA driver API 13.3**. Images must use a toolkit no newer than that ([CUDA guide, verified](https://modal.com/docs/guide/cuda)).
- **GPU strings:** `T4`, `L4`, `A10`, `L40S`, `RTX-PRO-6000`, … ([GPU guide](https://modal.com/docs/guide/gpu)).
- **Agent workflow:**
  - Install and authenticate: `uv tool install modal` (or `pip install modal`), then `modal token new`. That opens a browser sign-in and writes `~/.modal.toml`.
  - Run: `modal run tools/modal_gpu.py::test --gpu L4`. The script defines an image from `nvidia/cuda:13.0.2-devel-ubuntu24.04`, adds the local source tree, builds with CMake, then runs ctest and the benchmark. Output streams back to the Mac shell, and the exit code propagates.
  - Sandboxes (`modal.Sandbox.create(..., gpu=...)` plus `sb.exec(...)`) also work, up to 24 h ([docs](https://modal.com/docs/guide/sandboxes)). GPU Sandboxes require `runtime="gvisor"`, and CPU/memory cost 3× there, so prefer Functions.
- **Cold start (estimate, not measured):**
  - The first image build pulls ~4–5 GB of CUDA devel image, so it takes minutes.
  - After that the image is cached, and container plus GPU attach usually takes seconds to tens of seconds.
- **Gotchas:**
  - It runs inside gVisor. Nsight Compute hardware counters are probably unavailable (unverified), so plan on timer-based benchmarks.
  - Availability of L40S and RTX PRO 6000 can queue.
  - The free credit covers about 25 L4-hours, or about 20 h counting 4 CPU cores + 16 GiB RAM.
  - Whether a Starter account with no payment method hard-stops at $30 is unverified. Check this in the dashboard at sign-up.

### RunPod

- **Prices** ([pricing page, verified](https://www.runpod.io/pricing)):

  | GPU | Community | Secure |
  |---|---|---|
  | RTX 4090 | $0.34 | $0.89 |
  | **RTX 5090** | **$0.69** | $1.19 |
  | L4 | $0.44 | $0.59 |
  | L40S | $0.79 | $1.09 |
  | RTX PRO 6000 | $1.69 | $2.49 |

- **Cheaper Ada cards (tracker):** RTX 4000 Ada from ~$0.18 and RTX 4080 SUPER ~$0.28 ([computeprices](https://computeprices.com/providers/runpod/gpus/rtx4000ada)). RTX 5080 was not found.
- **Billing:** per second.
- **Storage:** volume disk $0.10/GB-mo while running and $0.20 while idle.
- **Windows:** no evidence of Windows pods. Pods are Linux containers.
- **Agent workflow:**
  - Auth: create an API key on the Credentials page (restricted scope is recommended), then run `runpodctl config --apiKey …` or `runpodctl doctor`.
  - Create a pod: `runpodctl pod create --gpu-id "NVIDIA GeForce RTX 5090" --cloud-type COMMUNITY --image nvidia/cuda:13.0.2-devel-ubuntu24.04 --min-cuda-version 13.0 --container-disk-in-gb 30`.
  - Then `runpodctl pod get <id>`, SSH in, and finish with `runpodctl pod delete <id>` ([CLI ref](https://docs.runpod.io/runpodctl/reference/runpodctl-create-pod)).
- **Gotchas:**
  - Community pods may have no public IP for SSH.
  - The account must hold at least one hour of credit to deploy ([billing docs](https://docs.runpod.io/get-started/billing-information)). A $10 minimum top-up is reported by third parties (unverified).
  - `ncu` hardware counters fail in unprivileged containers (`ERR_NVGPUCTRPERM`), per a [community report](https://www.answeroverflow.com/m/1225015304682471454).
  - Always pass `--min-cuda-version` so the host driver can run CUDA 13 binaries.

### Vast.ai

- **What it is:** a peer-to-peer marketplace; prices set by hosts (tracker figures).
  - RTX 4090: ~$0.20–0.48 ([Beam](https://www.beam.cloud/blog/nvidia-rtx-4090-pricing)).
  - RTX 5090: from ~$0.34 ([computeprices](https://computeprices.com/providers/vast/gpus/rtx5090)).
  - **RTX 5070 Ti: $0.13–0.24** ([Vast page](https://vast.ai/pricing/gpu/RTX-5070-TI)).
  - RTX 5080: from ~$0.17 ([computeprices](https://computeprices.com/providers/vast/gpus/rtx5080)).
  - RTX 5060 Ti: from ~$0.08 ([computeprices](https://computeprices.com/providers/vast/gpus/rtx5060ti)).
- **Billing:**
  - GPU is billed per second while running.
  - **Storage is billed while stopped**; only destroying the instance ends it.
  - Bandwidth is billed per byte at host-set rates ([billing](https://docs.vast.ai/billing)).
  - Requires prepaid credit.
- **Windows:** no. Docker instances, plus Linux VMs since Dec 2024 ([VM announcement](https://vast.ai/press-release/vast-ai-adds-virtual-machine-support-gpu-rental)).
- **Agent workflow** ([CLI docs](https://docs.vast.ai/cli/get-started)):
  - Install and auth: `pip install vastai`, then `vastai set api-key KEY`.
  - Find an offer: `vastai search offers 'gpu_name=RTX_5070_Ti num_gpus=1 verified=true rentable=true cuda_vers>=13.0' -o 'dph+'`.
  - Create: `vastai create instance <offer> --image nvidia/cuda:13.0.2-devel-ubuntu24.04 --disk 30 --ssh --direct`.
  - Use and tear down: `vastai ssh-url <id>`, `vastai copy`, `vastai destroy instance <id>`.
- **Gotchas:**
  - Host quality varies (PCIe width, thermals, noisy neighbours). Use `verified=true`, and record the host ID with every benchmark.
  - The `cuda_vers` filter field name is unverified; check it with `vastai search offers --help`.

### TensorDock

- **Prices:** marketplace KVM VMs. RTX 4090 from $0.35/hr on-demand ([TensorDock](https://marketplace.tensordock.com/gpu-4090.html)); older pages quote ~$0.20 spot.
- **Windows:** **yes, Windows 10 is a selectable OS**, with RDP needing ports 22 and 3389 forwarded and ≥80 GB disk ([RDP guide](https://docs.tensordock.com/virtual-machines/how-to-rdp-into-your-instance)).
- **API:** yes ([docs](https://docs.tensordock.com/llms-full.txt)).
- **Gotchas:**
  - Windows licensing cost, Windows GPU stock and billing granularity were not confirmed. Check in the console before renting.
- **Use:** the cheapest real Windows + RTX box, for an occasional Windows runtime smoke test.

### Lambda

- **Prices:** no consumer GPUs. A10 is listed variously at $0.69–1.29/hr ([costbench](https://www.costbench.com/compare/lambda-labs-vs-paperspace/), [computeprices](https://computeprices.com/providers/lambda/gpus/a10)). The rest of the range is H100/B200-class.
- **Windows:** Linux only.
- **Verdict:** not a fit. It has no Ada or Blackwell consumer parts, and its pricing is no better than Modal's.

### Paperspace / DigitalOcean GPU Droplets

- **Prices:** RTX 4000 Ada ~$0.76/GPU-h (secondary source); RTX 6000 Ada and L40S $1.57 ([DO pricing](https://www.digitalocean.com/pricing/gpu-droplets)).
- **Windows:** no verified Windows GPU option found.
- **Verdict:** not a fit. It costs more than RunPod or Vast.ai for similar parts.

### Google Colab

- **CLI:** an official CLI exists: [`googlecolab/google-colab-cli`](https://github.com/googlecolab/google-colab-cli), April 2026, actively maintained.
  - `colab run --gpu T4 script.py` provisions, runs and tears down.
  - `colab ssh` opens a shell.
  - It runs on Linux/macOS hosts and uses ADC/OAuth auth.
- **Free tier:** effectively a T4 only, with unpublished dynamic limits and 12 h max.
- **Restriction:** the [FAQ](https://research.google.com/colaboratory/faq.html) disallows "remote control such as SSH shells" on free runtimes without compute units.
- **Driver:** unknown; check with `nvidia-smi` before relying on CUDA 13.
- **Verdict:** fine for an occasional manual T4 sanity run. It is not a dependable agent CI target on the free tier.

### Kaggle

- **Free GPUs:** T4×2 or P100, ~30 GPU-h/week; figures are from secondary sources ([gpuperhour](https://gpuperhour.com/blog/free-cloud-gpus-and-credits)).
- **Automation:** possible via `kaggle kernels push`, but notebook-shaped and slow to start. A phone-verified account is required (unverified).
- **P100:** it is sm_60, which **CUDA 13 cannot target**.
- **Verdict:** not worth the friction.

### AWS, Azure, GCP Windows GPU VMs

- **AWS g6.xlarge** (L4, 4 vCPU, 16 GiB):
  - Linux on-demand $0.805/hr ([Vantage](https://instances.vantage.sh/aws/ec2/g6.xlarge)). Windows ~$0.99/hr, from an unlabeled aggregator column ([cloudpricecheck](https://cloudpricecheck-com.pages.dev/aws/ec2-pricing)), unverified.
  - Spot ~$0.35–0.62 (Linux) ([DoiT](https://compute.doit.com/spot/us-east-1/g6.xlarge)).
  - Windows needs the NVIDIA driver installed from AWS's driver bucket.
  - New accounts usually have **0 G-instance vCPU quota**, so a quota request is needed (common, unverified for this account).
  - Driven by `aws ec2 run-instances` plus SSM Run Command for PowerShell.
- **Azure:**
  - **Avoid `NV6ads_A10_v5`.** It is 1/6 of an A10 behind a vGPU (GRID) driver, and CUDA use is [questioned by Microsoft staff](https://learn.microsoft.com/en-gb/answers/questions/1659628/is-the-nvadsa10-v5-series-impossible-to-use-stable).
  - `NC4as_T4_v3` (full T4) Windows on-demand ~$0.71–0.75/hr; Linux spot ~$0.15–0.22 ([economize](https://www.economize.cloud/resources/azure/pricing/virtual-machine/nc4ast4v3/)).
- **GCP `g2-standard-4`** (L4):
  - Spot ~$0.17–0.43 ([gcloud-compute](https://gcloud-compute.com/g2-standard-4.html)).
  - Windows Server license adds ~$134/month-equivalent, about $0.18/hr. GPU quota is again 0 by default.
- **Verdict:** this is the only *managed* way to run the Windows build on an NVIDIA GPU. The GPUs are data-centre parts (L4/T4), though, and quota approval and Windows boot take time. Use it only if TensorDock or the user's own machine is unavailable.

### GitHub Actions

- **Standard runners** ([runner pricing, verified](https://docs.github.com/en/billing/reference/actions-runner-pricing)):
  - Windows 2-core is $0.010/min; Linux 2-core is $0.006/min.
  - **Free for public repos.**
  - For private repos, the Free plan includes **$12/month** of compute, about 1,200 Windows minutes ([billing](https://docs.github.com/en/billing/concepts/product-billing/github-actions)).
- **Runner images:**
  - `windows-latest` and `windows-2025` now map to the **VS 2026** image. **`windows-2022`** has VS 2022 17.14 ([runner-images](https://github.com/actions/runner-images)).
  - So pin `windows-2022` to match "MSVC 2022".
- **GPU runners:**
  - NVIDIA T4, 16 GB: `linux_4_core_gpu` $0.052/min ($3.12/hr) and `windows_4_core_gpu` **$0.102/min ($6.12/hr)**.
  - Larger runners need an org plan (Team/Enterprise, secondary source) and are **not free even for public repos**.
  - Verdict: expensive for a T4. The Windows GPU runner is the only zero-ops "Windows + NVIDIA" option, but at ~$6/hr.
- **Self-hosted runners:** free. The proposed $0.002/min fee was [postponed](https://www.theregister.com/2025/12/17/github_charge_dev_own_hardware/) and is not in current docs.
- **CUDA install action:** `Jimver/cuda-toolkit` v0.2.36 (2026-08-02) has Windows installers for CUDA 13.0.0 – **13.3.1**, but not 13.4. Checked via the GitHub API on `src/links/windows-links.ts`.

## 3. Can the Windows CUDA build be validated without a Windows GPU?

**Yes, for compilation. No, for runtime.**

- **Cross-compiling from macOS is impossible.** NVIDIA has shipped no macOS-host nvcc since CUDA 10.2, and Windows nvcc requires `cl.exe` as its host compiler.
- **Compiling needs no GPU or driver.** nvcc only needs the toolkit and MSVC. A `windows-2022` GitHub runner with CUDA installed compiles and links the whole project, including statically linked `cudart`. Tests that need a device must be labelled so CI can skip them (e.g., a ctest label `gpu`).
- **The Linux build can also run on the Mac.** NVIDIA publishes arm64 (SBSA) `nvidia/cuda` devel images, so nvcc runs in Docker on Apple Silicon. Device code is identical whatever the host CPU, so this catches nearly all `.cu` errors in seconds. It does not catch MSVC-specific host-code issues; the Windows CI job covers those.
- **This was tested here** with `docker run --platform linux/arm64 nvidia/cuda:13.0.2-devel-ubuntu24.04`. A C++20 kernel was built for `sm_75/86/89/120 + compute_120`; the result is in the Evidence section.
- **Runtime on Windows needs real hardware.** That means the user's machine, a TensorDock Windows VM, a Windows cloud VM, or the GitHub Windows GPU runner. Kernels behave the same on Linux and Windows with the same driver branch. The Windows-specific runtime risks are host-side: paths, the WDDM timeout (TDR) on long kernels, and driver version. An occasional Windows smoke run is enough, not one per commit.

## 4. Toolkit, driver and fatbin targets

### Toolkit facts

- **Current toolkit:** CUDA 13.4 Update 2 ([release notes](https://docs.nvidia.com/cuda/cuda-toolkit-release-notes/)).
- **Architectures:** CUDA 13.0 removed offline compilation for Maxwell, Pascal and Volta. 13.x supports **Turing (sm_75) through Blackwell** ([13.0 notes](https://docs.nvidia.com/cuda/archive/13.0.0/cuda-toolkit-release-notes/index.html)). sm_120 needs CUDA ≥ 12.8.
- **MSVC 2022 (193x):** supported in 13.0 through 13.4. VS 2026 (195x) is supported in at least 13.3 and 13.4 ([13.3 install guide](https://docs.nvidia.com/cuda/archive/13.3.1/cuda-installation-guide-microsoft-windows/index.html), [13.4](https://docs.nvidia.com/cuda/cuda-installation-guide-microsoft-windows/)). VS 2017 was dropped.
- **C++20:** listed as a supported dialect for MSVC 2022/2026.
- **Driver:**
  - The CUDA 13.x minor-version-compatibility floor is **driver ≥ 580**.
  - The Windows display driver is **no longer bundled** with the toolkit (Windows since 13.1). The user installs it from nvidia.com.
  - RTX 50 GPUs require a Blackwell-capable driver anyway.

### Recommendation

1. **Pin CUDA 13.0.x (13.0.2) everywhere:** the Mac Docker image, Linux CI, Windows CI (Jimver), Modal, RunPod and Vast.ai. 13.0 is the oldest 13.x, so the minimum driver is cleanly **R580**, with no minor-version-compatibility caveats. Modal's host (13.3) and `--min-cuda-version 13.0` filters cover it. Upgrade deliberately later. A 13.3 build would also run on R580 for SASS, but its PTX could not JIT on R580.
2. **Minimum driver on the user's machine:** NVIDIA **R580 or newer** (Windows 580.xx, Aug 2025). At startup, call `cudaDriverGetVersion` and refuse with a plain message if it is below `13000`.
3. **Ship targets:** `CMAKE_CUDA_ARCHITECTURES="75-real;86-real;89-real;120"`.

   | Entry | Emits | Covers |
   |---|---|---|
   | `75-real` | SASS sm_75 | RTX 20, GTX 16 |
   | `86-real` | SASS sm_86 | RTX 30 (also A10) |
   | `89-real` | SASS sm_89 | RTX 40 (also L4/L40S) |
   | `120` | SASS sm_120 + PTX compute_120 | RTX 50 (also RTX PRO 6000). The PTX lets future GPUs (sm_12x and later) JIT on a newer driver |

   - Within a major version SASS runs forward (8.6 → 8.9, 12.0 → 12.1), but per-arch SASS is tuned better, so keep both 86 and 89.
   - Don't use `sm_120a`/`f` arch-specific targets in shipped binaries; they are not forward compatible.
   - This build does **not** run on A100/H100 (sm_80/90). Add `80-real` only if data-centre runs are ever needed.
   - For dev iteration, build one arch: `-DCMAKE_CUDA_ARCHITECTURES=89`, or `native` on a GPU box.

## 5. Recommended setup (cheapest that satisfies a, b, c)

### What it costs per month (~20 GPU-hours)

| Use | Where | Hours | Rate | Cost |
|---|---|---|---|---|
| Correctness tests + Ada perf proxy | Modal L4 (4 CPU, 16 GiB) | 15 | ~$1.12/h incl. CPU/RAM | ~$17 → **covered by $30 free credit** |
| Upper Ada bandwidth class | Modal L40S | 2 | ~$2.27/h | ~$4.50 → covered (total ~$21.5 < $30) |
| Blackwell sm_120 perf | Vast.ai RTX 5070 Ti | 3 | ~$0.13–0.24/h + storage/bandwidth | **~$1** |
| Compile checks | Mac Docker + GitHub `ubuntu-24.04` / `windows-2022` | — | free for a public repo; private repo ~$0.01/min Windows against $12/month included | **$0** (private: $0 until ~1,200 Windows min/month) |
| Windows runtime smoke | User's machine (free) or TensorDock Windows RTX 4090 | ~1 | ~$0.35/h | **$0–0.50** |
| **Total** | | ~20 | | **≈ $0–5/month cash** (Vast.ai needs a prepaid top-up, e.g. $10, which lasts months) |

Fallback if Vast.ai stock or host quality is poor: RunPod RTX 5090 Community at $0.69/h, about $2 for 3 h. Scale its results down by bandwidth (1792 → 896 GB/s) and label them as projections.

### Step-by-step (agent does everything after each sign-in)

**Spend gate:** the owner approves any card or top-up. Nothing below costs money until a top-up is made, except usage beyond Modal's free credit.

1. **Local compile check (no sign-in).**
   - Add a script, e.g. `scripts/cuda-compile-check.sh`. It runs `docker run --rm --platform linux/arm64 -v "$PWD":/w -w /w nvidia/cuda:13.0.2-devel-ubuntu24.04` with `cmake -G Ninja -DCMAKE_CUDA_ARCHITECTURES="75-real;86-real;89-real;120"` and a build.
   - Agents run it after every CUDA change. Docker is already installed on this Mac (server 29.8.0, arm64).
2. **GitHub CI.** The owner creates or approves the GitHub repo (public is free; `gh` is already signed in as `ezzy1630`). The agent adds a workflow with two jobs, both compile-only and both running non-GPU unit tests:
   - `ubuntu-24.04`: `Jimver/cuda-toolkit` with `cuda: '13.0.2'`, `method: network`, sub-packages `nvcc`, `cudart`. Then CMake build and `ctest -LE gpu`.
   - `windows-2022`: `ilammy/msvc-dev-cmd` (or the VS dev shell), then `Jimver/cuda-toolkit` with `cuda: '13.0.2'`, `method: network`, sub-packages `["nvcc","cudart","visual_studio_integration"]`, `use-github-cache: true`. Then a CMake (Ninja) Release build and `ctest -LE gpu`.
3. **Modal (sign-in 1).**
   - The agent runs `uv tool install modal && modal token new`. The owner completes the browser sign-in with GitHub or Google and **does not add a card**.
   - Then check in the dashboard whether usage hard-stops at the $30 credit.
   - The agent adds `tools/modal_gpu.py`. It defines an image from `nvidia/cuda:13.0.2-devel-ubuntu24.04` with cmake/ninja, uses `add_local_dir` for the source, and has `test` and `bench` functions taking a GPU type (`L4` default; `L40S`, `T4` optional via `with_options`).
   - Agents run `modal run tools/modal_gpu.py::test` per meaningful change, and `::bench` on demand. Results come back as JSON on stdout: GPU name, driver, measured STREAM BW, cells/s, %BW, commit SHA.
4. **Vast.ai (sign-in 2, plus spend approval).**
   - The owner signs up and approves a small prepaid top-up. The agent runs `pip install vastai && vastai set api-key …` and registers a local SSH public key.
   - Each benchmark session follows the same script: `search offers 'gpu_name=RTX_5070_Ti verified=true rentable=true …' -o dph+`, then `create instance` with the CUDA 13.0.2 image, then `copy` the source, build, run the bench, `copy` the results back, and `destroy instance`.
   - **Always destroy, never stop**: stopped instances keep billing storage.
5. **Windows runtime smoke (per release, not per commit).**
   - Preferred: with the user's consent, run the release `.exe` on their own RTX machine. Later, optionally register it as a GitHub self-hosted runner (free, exact hardware) with the label `windows-rtx`.
   - Alternative: rent a TensorDock Windows 10 RTX 4090 VM for an hour.
6. **Benchmark hygiene.**
   - Every result records the provider, host ID, GPU name, driver version, CUDA version, commit, measured STREAM bandwidth and clocks (`nvidia-smi -q -d CLOCK`).
   - Laptop numbers are labelled as projections.

## Evidence and limits

- **Tested on this Mac:** arm64 and amd64 variants exist for `nvidia/cuda:13.0.2-devel-ubuntu24.04` and `13.3.1-devel-ubuntu24.04` (`docker manifest inspect`).
- **Docker compile probe, run 2026-10-29:** `docker run --platform linux/arm64 nvidia/cuda:13.0.2-devel-ubuntu24.04` with `nvcc -std=c++20` and gencodes for 75/86/89/120 plus `compute_120` PTX.
  - nvcc reported `release 13.0, V13.0.88`.
  - `cuobjdump` listed `sm_75`, `sm_86`, `sm_89` and `sm_120` cubins plus `sm_120.ptx`; exit 0.
  - First run took 1:37 including the image pull. The container warns "NVIDIA Driver was not detected", which is expected because compiling needs none.
- **Read via the GitHub API on 2026-10-29:**
  - Jimver's supported Windows CUDA versions.
  - The VS versions on GitHub's runner images (`windows-2022`: VS 2022 17.14; `windows-latest` → `windows-2025-vs2026`).
- **Read from vendor pages on 2026-10-29:** prices for Modal and RunPod, GitHub runner prices, Modal's host driver version, and the CUDA 13.0/13.3/13.4 release notes and install guides.
- **Unverified:**
  - Marketplace prices (Vast.ai, TensorDock), which move daily.
  - The Windows rates for AWS, Azure and GCP.
  - The free-tier limits for Colab and Kaggle.
  - Whether Nsight Compute counters are available on Modal.
  - The plan requirement for GitHub GPU runners.
  - All cold-start times, which are estimates and were not measured.
