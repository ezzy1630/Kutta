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
  --approximate SETTING[,SETTING...]
                               Explicitly permit each named approximation.
                               Currently only wfgc is allowed; no wildcard.

Outputs
  result.vtu      Cell fields u (velocity vector), p (gauge pressure), k,
                  omega, mu_t; open with ParaView. No boundary-face fields:
                  the current VTU writer supports cell fields only.
  forces.csv      Pressure/viscous/total wall forces in N, report projections
                  and coefficients for resolved force/drag/lift definitions.
  residuals.csv   Fluent-scaled residuals at solver monitor intervals, plus
                  average wall seconds per iteration for each interval.
  summary.json    SHA256 input hashes, version/commit, device, cell/iteration
                  counts, setup/solve seconds, stopping/convergence status,
                  approximate flag and the approximations actually used.

Limitations
This is an early solver, not yet validated as a Fluent replacement. By default,
it accepts only settings the strict adapter can reproduce; unsupported settings
are listed and exit nonzero. Do not change your Fluent case merely to hide a
rejection.

Optional approximation (changes the discrete equations)
  .\kutta.exe solve "C:\cases\car.cas.h5" --backend cuda --approximate wfgc --out "C:\cases\kutta-approximate"

wfgc permits recon/cell-lsf? = #t (warped-face gradient correction on) to run
with WFGC off. Every other unsupported setting still fails, and an absent or
unparsed WFGC setting still fails. Unknown, empty and duplicate names fail.
The case file is never changed. Without --approximate, WFGC on still fails
and the error names the opt-in. Permission has no effect when WFGC is off.

Effect documented in src/kutta/fv/VERIFICATION.md, "V6 initialization comparison":
the CPU hybrid mixing-elbow diagnostic with WFGC off differed from Fluent by
3.366748% velocity L2, 2.665514% pressure L2, 14.722368% k L2, 2.805466% omega
L2 and 1.973008% pressure drop. These are combined solver discrepancies, not
isolated WFGC errors or bounds for your case. This is an approximate diagnostic,
not a faithful reproduction or an accepted same-settings accuracy result.

When used, APPROXIMATE RESULT is printed at start and end. summary.json records
approximate=true and each setting, Fluent value, Kutta action, measured effect
and documentation path. Both CSV files begin with a '# APPROXIMATE RESULT'
comment; skip lines starting with # when reading CSV. result.vtu carries the
same warning in an XML header comment. Keep these markers when sharing results.
With no approximation used, CSV/VTU headers are unchanged and the summary has
approximate=false and an empty approximations list.

The adapter does not yet expose convergence criteria: the CLI stops
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
