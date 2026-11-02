Evidence for the native tet/prism convergence investigation.

`baseline_reproduction.json` compares every non-timing history column with the
previous lane's retained candidate histories. Metal matches through 5,000
iterations; the instrumented CPU prefix matches through 1,500.

`cpu_cycle_cells.csv` contains 100 consecutive step snapshots for five selected
input cells. IDs are zero-based input IDs before RCM. `cpu_simplec_trace.json`
summarizes all 130,331 cells over iterations 1,000 through 1,099.
`cpu_residual_localization.json` records concentration of continuity and the
iteration-1,500 momentum residuals. The pressure-force column is the local
Gauss pressure source, divided by cell volume. It is not a body force integral.

Field timing matters. Step velocity, pressure, k and omega are the updated
state. Gradients, limiters and assembled sources belong to the preceding
assembly. `continuity` is the provisional predictor imbalance before pressure
correction, matching the solver's reported raw-continuity contract. Monitor
snapshots refresh the momentum equation assembly at the current state. The
trace only downloads fields; it does not call an extra monitor or change
limiter state. Full binary dumps remain under ignored
`runs/tet_prism_convergence/trace_cpu_simplec`.

The binary dump has 21 float32 arrays, each with one value per RCM cell, in
this order: u, v, w, pressure, k, omega, u/v/w limiter, pressure-gradient x/y/z,
momentum-source x/y/z, provisional continuity, physical momentum diagonal,
pressure correction, and pressure-source x/y/z. Monitor dumps append three
absolute momentum-residual arrays. `analyze_trace.py` consumes this layout.
A strongest FFT bin spanning all 100 iterations means unresolved drift or a
longer period; it does not establish a 100-iteration oscillation.

The linear and box probes test pressure reconstruction independently of
native drag. For the linear probe, mode 0 is the original fit, mode 1 moves
only amplified symmetry fallback rows to the normal foot, and mode 2 gives
all symmetry faces normal-foot rows. The stronger exact symmetry-constraint
control is separately named. The normal-foot correction is the only retained
solver experiment, and remains opt-in.

`option_regression_red.log` was captured after adding the option to the type
but before using it in the fit. `option_regression_green.log` records the
first passing CPU/Metal flow regression. Later verification also checks that
the option rejects the unsupported Green-Gauss combination. No bound or
iteration count was changed after the red result.

Rejected prototype runs are recorded separately from retained-build runs.
Their embedded build stamp names the base commit, not the uncommitted
experimental source. Raw histories contain wall-clock columns for audit;
no performance conclusion uses those timings. Failed runs do not supply an
accepted drag coefficient or a physical time average.

`poly_mesh.json` records convex tet-pair agglomeration, exact coordinate and
boundary-role readback, cell-volume preservation and refusal to overwrite an
existing output. This is a generated polyhedral variant of the native mesh,
not real Fluent Meshing output. The large mesh remains ignored.

Slow XML and restricted fixtures are not committed. The serial slow runner
changes only the process boundary around each selected case. It retains all
37 cases, backends, subcases, solver budgets and assertions. Per-case resource
and assertion receipts bind the merged XML to the binary. Original quick
cases outside the selected slow suite are not included in the merged XML.
The shared runs root lacks two generated wing files; the supplemental
conversion case uses locally regenerated files with the exact ledger hashes.

`final_source.json` binds the clean retained source and frozen binaries at
`998a666` to the final runs. `historical_b27_source.json` preserves the earlier
symmetry-only implementation's provenance. Rejected prototypes retain their
base build stamps in their logs; those stamps alone do not identify the scratch
source. The combined-fallback and final retained histories match all 100 samples
exactly (`retained_reproduction.json`).

`boundary_regression_red.json` and `boundary_regression_green.json` preserve the
actual expressions and backend/boundary contexts without committing test XML.
Their receipts include source and binary hashes. The older `option_regression_*`
logs cover the initial symmetry-only correction. The later red result passes
symmetry and fails moving walls; the final option passes both.

`constrained_turbulence_localization.json` describes the rejected exact-symmetry
Coupled control, not the final option. Its 250-iteration scalar snapshots locate
k and omega residuals and record source/sink coefficients; they cannot resolve
an oscillation frequency or isolate a causal turbulence term. `strain_affine.log`
and `strain_control_source.json` bind the separate rejected strain control.

The final inventory contains 35 completed benchmark records, including the
1,500-iteration CPU diagnostic and runs that stopped with a nonfinite failure.
Tail-window coverage and sample counts are explicit; failed runs retain null
means and null accepted force values. These records do not turn an iteration
cap into convergence.

Three later controls are `projected_strain_metal_simplec`,
`diff_velocity_normal_foot_metal_simplec` and
`diff_turbulence_normal_foot_metal_simplec`. The first is a host-side strain
trace prototype; the second uses only retained library options; the third
adds smooth SST scalar limiting with the same existing coefficient. None meets
the full gate, and neither prototype is retained. Their source receipts and
launch records distinguish them from their embedded base commit stamps.
The launch records also show pauses/resumptions of this lane's iteration-limited
Tyrrell test while these controls ran serially. Pauses affect elapsed time only;
completed slow-case assertion values are compared with the untouched baseline.

`interruption_recovery.json` records the interrupted slow case 23 and the
runner's resume guards. The first 22 completed case XML files remain immutable.
The interrupted case's partial XML is preserved only under ignored runs and
is excluded from acceptance. Resume starts that entire case from its initial
state with its original budgets; it does not continue a solver state or count
partial grid results as a completed case. The final comparison must include
all 37 complete case receipts.

The pressure-correction-only control and two combined-fit controls are also
complete. `correction_fit_probe.json` isolates the tiny-case recovery with
physical pressure unchanged. All three full native controls fail the original
5,000-iteration gate. Their launch records and `correction_foot_source.json`
bind the numerical prototype; it is not retained. The final slow comparison
contains all 37 complete cases, with zero changed assertion values against the
fresh 22e65c7 baseline and exactly the prior 31 changes against 8866cc0.
