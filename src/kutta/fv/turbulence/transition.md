# Fluent Intermittency Transition Model

Kutta implements Fluent's **Intermittency Transition Model** (one-equation
gamma, Menter et al. 2015, Fluent's `gamma-transport-eqn`) on top of SST.
Fluent's Transition SST (gamma-Re_theta, four equations) and its algebraic
transition model are rejected by the case adapter with a reason.

## Why this model

Fluent 2026 R1 prints the intermittency model completely: the gamma transport
equation, every onset and pressure-gradient function, and every default
constant (Theory Guide 4.8.2), plus the k-equation coupling it shares with
Transition SST (4.8.3, which points to 4.7.2.2, Eqs. 4-205..207). The Transition
SST page instead gives `Re_theta_t = f(Tu, lambda)`, `F_length = f(...)` and
`Re_theta_c = f(...)` (Eq. 4-201). It refers to Langtry and Menter for the
correlations and says that Re_theta_t "has been modified", without printing the
modification. Those correlations cannot be reproduced exactly from the
documentation.

## Equations (Fluent TG 4.8.2)

```
d(rho g)/dt + d(rho U_j g)/dx_j = P_g - E_g + div[(mu + mu_t/sigma_g) grad g]
P_g = F_length rho S g (1 - g) F_onset,              F_length = 100
E_g = c_a2 rho Omega g F_turb (c_e2 g - 1),           c_a2 = .06, c_e2 = 50
F_onset1 = Re_V / (2.2 Re_theta_c), F_onset2 = min(F_onset1, 2)
F_onset3 = max(1 - (R_T/3.5)^3, 0), F_onset = max(F_onset2 - F_onset3, 0)
F_turb = exp(-(R_T/2)^4), R_T = rho k/(mu omega), Re_V = rho d^2 S/mu
Re_theta_c = C_TU1 + C_TU2 exp(-C_TU3 Tu_L F_PG(lambda_thetaL))
Tu_L = min(100 sqrt(2k/3)/(omega d), 100)
lambda_thetaL = -0.1111 (dV/dy) d^2/nu + 0.1875, bounded to [-10, 10]  (4-216, 4-217)
F_PG = min(1 + C_PG1 lambda, C_PG1_lim)                            lambda >= 0
F_PG = min(1 + C_PG2 lambda + C_PG3 min(lambda + 1, 0), C_PG2_lim) lambda < 0
```

Defaults: C_TU1 = 100, C_TU2 = 1000, C_TU3 = 1, C_PG1 = 1, C_PG2 = -0.5,
C_PG3 = 0, C_PG1_lim = 1.5, C_PG2_lim = 3, C_SEP = 1, sigma_gamma = 1. Fluent
rescales lambda by 14.68 relative to the paper, so these values are the paper's
C_PG1 = 14.68 and C_PG2 = -7.34.

k-equation coupling (TG 4-205..207):
`G*_k = gamma G~_k` and `Y*_k = min(max(gamma, 0.1), 1) Y_k`. The omega
equation is unchanged. G~_k is SST's production after the case's production
option, its limiter and the wall treatment, so Kutta follows the case's
production setting. The paper uses Kato-Launder; Fluent documents the original
production term.

## Terms Fluent names but does not print

These come from Menter et al. 2015 as restated on NASA TMR's
[SST-2003 gamma page](https://tmbwg.github.io/turbmodels/menter_gamma_3eqn.html)
(formerly turbmodels.larc.nasa.gov):

- `P_k^lim = 5 C_k max(gamma - 0.2, 0)(1 - gamma) F_on^lim max(3 C_SEP mu - mu_t, 0) S Omega`,
  with `F_on^lim = min(max(Re_V/(2.2*1100) - 1, 0), 3)`. It is added to the k
  production. Fluent names C_SEP in the user interface and the case stores
  `trans-one-eqn-csep`.
- `F1 = max(F1_SST, F3)`, `F3 = exp(-(R_y/120)^8)`, `R_y = rho d sqrt(k)/mu`.
  The k and omega diffusivities that depend on F1 are recomputed. mu_t depends
  only on F2 and is unchanged.
- `dV/dy = n . grad(n . U)`, with n the unit vector from the nearest wall point
  (`wall_direction()` in `device_mesh.hpp`). This is TMR's general-geometry form
  and Fluent's `trans-one-eqn-dudn-method 1`; other methods are rejected.
- Boundary conditions: zero gamma flux at walls, and the zone's `intermit` value
  on inflow faces. Fluent cases store `(intermit (constant . 1) ...)` on inlet
  and outlet zones.
- Robustness floor `max(F_PG, 0)`. It cannot bind with C_PG3 = 0.

## Inferences

- `trans-one-eqn-tke-term-coef` is taken to be C_k of P_k^lim. The mapping is
  inferred from the name, and its default 1 matches the paper.
- The rpvar names (`rp-trans-one-eqn?`, `trans-one-eqn-*`, `intermit/*`,
  `trans-intermit/min-limit` and `max-limit`) and their defaults were read from
  validation/fluent_examples mixing-elbow cases, in which the model is
  **inactive**. No case with the model active was available, so the adapter's
  reading of an active case is untested against Fluent output.
- Fluent clamps gamma to `[trans-intermit/min-limit, max-limit]` (1e-7, 1) after
  each solve, as for its other limit rpvars.
- The source linearization is Kutta's choice: `P - E = (a+b) g - (a + c_e2 b) g^2`,
  with `a = F_length rho S F_onset` and `b = c_a2 rho Omega F_turb`. The positive
  part is explicit and the quadratic part becomes an implicit sink. Both
  coefficients are nonnegative, and the steady solution is unchanged.

## Not implemented (rejected with a reason)

- Crossflow transition (`trans-one-eqn-crossflow?`), intermittency smoothing,
  roughness correlations, and dudn methods other than 1.
- The model with any viscous model other than SST.

## Code

| Part | File |
| --- | --- |
| Kernels (gamma sources, diffusion, boundary, clamp, F1/F3, k coupling) | `src/kutta/kernels/transition.kern.h` |
| Gamma transport | `intermittency.hpp/.cpp` |
| k-equation coupling | `SstOptions::intermittency` in `sst.hpp/.cpp` |
| Fluent case parsing | `src/kutta/case/fluent_case.cpp` (`transition_model`) |
| Admission, constants, numerics, zone values | `src/kutta/setup/fluent.cpp` |
| Oracle and parity tests | `tests/turbulence/test_transition.cpp` |
| Validation (T3A/T3B/T3A-, S1223, TMR plate) | `validation/results/transition/REPORT.md` |
