# Hybrid initialization

Implemented 2026-10-30 for the strict steady, constant-density SST adapter.
Initialization changes the starting state and convergence path, not any steady
equation, boundary condition, material or turbulence coefficient.

Sources, read 2026-10-30:

- [Fluent User's Guide 37.11](https://ansyshelp.ansys.com/public/Views/Secured/corp/v252/en/flu_ug/flu_ug_sec_init_hybrid.html):
  iteration count, scalar explicit relaxation, reference frame, pressure from
  inlets, external-aero linear potential, constant speed and averaged turbulence.
- [Fluent Theory Guide 23.8](https://ansyshelp.ansys.com/public/Views/Secured/corp/v252/en/flu_th/flu_th_sec_hybrid_initialization.html):
  velocity potential, wall/inlet/outlet conditions, pressure recipe and turbulence.

## Publicly specified behavior

Solve `laplacian(phi)=0`, with velocity `grad(phi)`. Walls and symmetry have zero
normal velocity; velocity inlets supply their specified normal velocity; pressure
outlets have zero potential. Tangential no-slip is imposed by the subsequent
momentum solve, not by this inviscid auxiliary problem.

Pressure has a separate Laplace problem when pressure information is available
at an inlet and an outlet. A velocity inlet supplies pressure only when the
specified-inlet-pressure option is enabled. Otherwise initialize pressure to
an average of the available boundary pressures. Walls have zero normal pressure
gradient. The adapter still rejects unsupported physical boundaries, moving
frames, porous media, FMG and active equation settings.

The initialization plan is a `variant<CellFields, HybridInitialization>`.
`setup::initial_fields()` parses it before choosing a Device;
`setup::initialize(solver, plan)` executes the same path in the CLI and V6.
Saved data continues to use `CellFields` in original source-cell order.

## Stored controls and explicit implementation choices

The public guides do **not** specify a reproducible Fluent implementation.
The following choices are Kutta's; this is not a claim of bitwise Fluent parity.

| Incomplete published detail | Kutta choice |
| --- | --- |
| Laplacian discretization, gradient reconstruction and AMG operations | Existing deterministic FV diffusion, independent LSQ boundary observations, stored-float M-matrix assembly and Kutta AMG-PCG. No host solve or atomics. Neumann observations are at the normal foot, not the displaced face centroid. |
| Meaning of an initialization iteration, residual normalization and linear tolerance | `hyb-init/niter` outer non-orthogonal/deferred-correction passes. Each linear solve uses the corresponding stored `hyb-init/residual-tolerance` as a relative L2 target, absolute target `1e-12`, maximum 300 Krylov iterations. Fluent's AMG cycle/method parameters are replaced by Kutta AMG, as reported by the CLI. This is initialization only. |
| Two scalar relaxation controls | Use the two entries of `hyb-init/relax` for implicit relaxation and `hyb-init/explicit-relax` for explicit blending with the preceding auxiliary iterate. Invalid/missing controls are rejected. |
| Auxiliary initial guesses and external-aero option's precise patch | Parse `hyb-init/patch`, `hyb-init/patch-ads?`, and `hyb-init/solve?`. An enabled linear patch uses the area-averaged inlet velocity, with potential origin at the area-averaged outlet centroid. Otherwise start from the stored constant scalar guess. A disabled scalar solve retains its initial guess; velocity still derives from that potential and the prescribed normal observations. |
| Float accuracy in a large domain | Solve a perturbation about the analytically known linear inlet-velocity potential, including its exact boundary fluxes and outlet values. Add that known velocity analytically after differentiation. This avoids subtracting nearly equal large float potentials across thin wall cells. The representation does not change the Laplace problem. |
| Mapping potential to conservative face mass flux | Retain the auxiliary diffusion operator's face flux, including internal and outlet non-orthogonal terms, for the first transport assembly. Differentiating a cell potential alone is not discretely solenoidal. Subsequent iterations use the ordinary steady flux operator. |
| Averaging measures and multiple inlet turbulence specifications | Area-average the inlet **recipes**, then evaluate them using local potential speed. Direct k/omega values remain constants; intensity/length gives k proportional to speed squared and omega proportional to speed; intensity/viscosity ratio gives both proportional to speed squared. With `hyb-init/turb-domain-averaged?`, volume-average the resulting k/omega fields; otherwise retain them locally. Floors are applied before averaging. |
| No inlet turbulence data | Use stored `hyb-init/turb-intens` and `hyb-init/turb-vr` as the local intensity/viscosity-ratio recipe. Closed stationary domains have zero potential velocity and the configured turbulence floors. |
| Constant speed's value and stagnation points | `hyb-init/constant-vel-vector?` selects the documented constant-magnitude option. Use the area-averaged inlet speed and the potential direction; leave a zero direction zero. Rescaling invalidates the potential's conservative flux, so use the ordinary flux reconstruction in that mode. |
| Pressure averaging measure | Area-average all eligible pressure boundaries. `hyb-init/pressure-from-bc?` additionally admits stored inlet `p`; otherwise only pressure outlets contribute. With none, use the stored pressure scalar guess. |
| Absolute versus relative reference frame | Equivalent in the adapter's supported stationary cell zones. Rotation/compressibility/species/multiphase recipes are outside that adapter's accepted physics. |

The real elbow stores `hyb-init? #t`, ten iterations, unit explicit and implicit
relaxation, `1e-6` scalar tolerances, averaged turbulence, and both scalar solves
enabled. Its standard defaults are intentionally not read for a hybrid start.
The old standard proxy remains available only as the explicit V6 diagnostic
`KUTTA_V6_STANDARD_PROXY=1`.

There is no separate "explicit averaging" switch in the inspected case or the
cited hybrid UI documentation. The documented explicit controls are scalar
under-relaxation; turbulence averaging is the separate boolean above.

The measured CPU/Metal initialization comparison and remaining original V6
accuracy/settings failures are recorded in `src/kutta/fv/VERIFICATION.md`.
