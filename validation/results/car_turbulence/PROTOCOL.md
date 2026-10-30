# V7 car turbulence protocol, frozen before measurements

SA reuses the V5 grids, dimensionalization, BC mapping, uniform start,
SOU momentum / FOU turbulence, 3,000 step / 600 s per-grid budget, CPU and Metal,
and the 100-step convergence sampling interval. No reference field seeding.
The reference is NASA standard SA including ft2; Fluent's documented no-ft2
variant is tested independently at the constitutive level. Smooth resolved
walls use nu-tilde=0 and molecular stress; this is not Fluent automatic walls.

| Gate | Pass condition |
|---|---|
| SA-1 equations | independent double equations, relative error <=3e-5 |
| SA-2 backend parity | all model fields CPU vs Metal, relative error <=3e-5; deterministic repeats |
| SA-3 BCs | direct nu-tilde, exact damped viscosity-ratio inversion, zero wall, inlet/backflow and outflow |
| SA-4 plate | finest Cf error <=5%, Cd <=7% vs finest CFL3D |
| SA-5 bump | each finest Cf error <=10%, each Cp absolute error <=.03, Cd error <=10% |
| SA-6 trend | >=3 levels, finest maximum reference-error norm < coarsest |
| SA-7 iterative | every level u/v/nu-tilde/provisional-continuity <=1e-5, load change <=.1% over 100 steps |
| SA-8 Fluent automatic wall | full calibrated blend plus independent external evidence; unavailable blend must throw |
| KE-1 realizable equations | independent double or published code comparison, variable Cmu and epsilon destruction |
| KE-2 walls | scalable log law including 11.225 floor; EWT two-layer plus enhanced wall functions |
| KE-3 flow | three-grid fully developed channel with MKM actual Re_tau=587.19: max U error <=10% for y+>=30, wall shear error <=2%, all equation defects <=1e-5 and last-100 U change <=1e-4; 10000 step / 180 s per-grid budget |
| Release | SA-1..8 and KE-1..3 where applicable, case roundtrip and strict unknown rejection, quick tier green |

Accuracy checks remain asserted for unconverged runs, but cannot establish
acceptance. These engineering allowances follow V5, not TMR's much smaller GCI.
No Fluent executable or same-mesh reference is available. Standard/RNG are
conditional on SA and realizable acceptance, and must remain disabled otherwise.
