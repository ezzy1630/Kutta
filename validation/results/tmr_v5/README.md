# V5 evidence

These are diagnostic artifacts, not accepted V5 solutions. The complete report
and reproduction commands are in `src/kutta/fv/VERIFICATION.md`, V5.

`measurements.csv` is generated from the final CPU/Metal slow-suite output with
`tools/report_tmr.py --test-log validation/results/tmr_v5/slow_v5.log`.
The raw `slow_v5.log` is retained beside the CSV. Runtime/failure reporting is
commit `1f12a3c`; the tested binary SHA-256 is
`c318750b118ee3e9ae009395654b7b633a5dce162b2db0523fee04448a633ef2`. It keeps failed
solutions as NaN and records their failure messages. `flatplate_profiles.csv`
uses the 137x97 CPU cell export, interpolated first along each original grid
row to x=0.970084071, then in y; the reference is the finest CFL3D profile.
These profiles compare an incompressible solution with M=0.2 reference data.

`bump_pressure_iteration18.coo.gz` is a captured, **unmodified** device float
pressure-correction matrix for the 89x41 bump (3,520 extruded cells). It was
captured immediately before `pressure_solver->refresh`, after 18 completed
outer steps, on the serial CPU backend. Configuration: U=rho=1, nu=1/3e6,
SST defaults, uniform initial fields k=2.25e-7 / omega=75, dt=.01,
momentum relaxation=.5, pressure relaxation=.3, span=.01, SOU momentum / FOU
turbulence. Numerics are commits `5ff0bf9` and `a1d6d7e`; setup is the
`validation/cases/tmr.hpp` implementation subsequently committed in `a8abab5`.
The temporary capture downloaded the assembled CSR values and printed each
stored float with 17 digits; its instrumentation has been removed. The first
line is cell count and completed iteration, followed by row/column/value
triples in RCM cell order. Gzip mtime is zero.

Compressed SHA-256:
`684171e2c30a9d05b55fb123e7f35f61b6738039c7286433cde2145bd58395a6`.

Recheck with `python3 tools/analyze_tmr_pressure.py` (NumPy required). Symmetry
error is exactly zero, but double-precision eigenanalysis gives a smallest
eigenvalue of approximately -1.10101827e-11 and Cholesky fails. Diagonals range
from 1.24645387463e-14 to .00592642463744. There are 1,712 negative row sums.
This demonstrates loss of positive definiteness in the stored fine operator;
it does not by itself isolate every later divergence in the finer bump grids.


`flatplate_convergence.svg` and `.png` are generated with
`python3 tools/plot_tmr.py` (Matplotlib required). Kutta points are explicitly
labelled unconverged. The figure includes same-grid CFL3D/FUN3D and
incompressible FUN3D trends, not only a single favourable finest-grid number.
