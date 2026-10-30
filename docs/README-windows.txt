Kutta for Windows (early development)
=====================================
Requirements
- Windows 10/11, 64-bit x86.
- NVIDIA RTX 20-series or newer; NVIDIA driver R580 or newer.
- Microsoft Visual C++ 2015-2022 x64 Redistributable:
  https://aka.ms/vs/17/release/vc_redist.x64.exe
- No CUDA toolkit, HDF5 installation, Python, or Fluent license is needed to run.

Extract the entire zip before running. Keep the supplied DLLs beside kutta.exe.
Inputs and results stay on your computer. Nothing is uploaded.

Examples (PowerShell, from the extracted directory)
  .\kutta.exe inspect "C:\cases\car.cas.h5"
  .\kutta.exe solve "C:\cases\car.cas.h5" --backend cuda --iterations 2000 --out "C:\cases\kutta-result"
  .\kutta.exe solve "C:\cases\car.cas.h5" --data "C:\cases\car.dat.h5" --backend cuda --iterations 2000 --out "C:\cases\kutta-result"

Options
  --backend cuda|metal|cpu|auto  Default auto: CUDA, Metal, then CPU.
  --iterations N               Positive iteration cap; default 1000.
  --out DIR                    Default kutta-result, relative to current folder.
  --data FILE.dat.h5           Seed fields from a matching Fluent data file.

Outputs
  result.vtu      Cell fields u (velocity vector), p (gauge pressure), k,
                  omega, mu_t; open with ParaView. No boundary-face fields:
                  the current VTU writer supports cell fields only.
  forces.csv      Pressure/viscous/total wall forces in N, report projections
                  and coefficients for resolved force/drag/lift definitions.
  residuals.csv   Fluent-scaled residuals at solver monitor intervals, plus
                  average wall seconds per iteration for each interval.
  summary.json    SHA256 input hashes, version/commit, device, cell/iteration
                  counts, setup/solve seconds, stopping/convergence status.

Limitations
This is an early solver, not yet validated as a Fluent replacement. It accepts
only the settings the strict Fluent adapter can reproduce; unsupported settings
are listed and exit nonzero. Do not change your Fluent case merely to hide a
rejection. The adapter does not yet expose convergence criteria: the CLI stops
at --iterations and records converged=false. Peak device memory is unavailable
and recorded as null. A successful exit means the requested run and outputs
completed; it does not establish convergence or Fluent force agreement.

Use a new output directory per run. Successful runs replace the four named
output files. A failed run may leave partial files; summary.json is written last.

Bundled dependencies
HDF5 1.14.6 and zlib 1.3.1, from vcpkg baseline
84bab45d415d22042bd0b9081aea57f362da3f35 (2025.12.12).
CUDA runtime is linked statically. Third-party license notices are in licenses/.
The executable source commit is recorded in build-info.txt and summary.json.
