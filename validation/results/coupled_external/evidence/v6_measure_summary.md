
### Fluent residual criteria (kutta solve)

| case | backend | Coupled before | Coupled after | SIMPLEC |
|---|---|---|---|---|
| elbow | cpu | 160 / 15.1–15.3 s | 160 / 15.5–15.7 s | 240 / 26.2–26.4 s |
| elbow | metal | 160 / 2.2–2.3 s | 160 / 2.2–2.3 s | 260 / 3.8–3.9 s |
| nozzle | cpu | 60 / 7.4–7.4 s | 60 / 7.7–7.7 s | 80 / 10.2–10.3 s |
| nozzle | metal | 60 / 1.2–1.3 s | 60 / 1.2–1.3 s | 80 / 1.2–1.3 s |

### V6 criterion (bench_convergence)

| case | backend | Coupled before | Coupled after | SIMPLEC |
|---|---|---|---|---|
| elbow | cpu | 230 / 24.0–24.1 s | 220 / 23.3–23.5 s | 500 / 62.2–62.4 s |
| elbow | metal | 230 / 3.9–4.1 s | 220 / 3.7–3.8 s | 500 / 9.7–10.0 s |
| nozzle | cpu | 70 / 9.2–9.3 s | 70 / 9.4–9.5 s | 150 / 21.2–21.3 s |
| nozzle | metal | 70 / 1.6–1.6 s | 70 / 1.6–1.6 s | 150 / 2.9–3.1 s |

### Coupled solver bytes per cell (CoupledSolver::bytes)

| case | before | after |
|---|---|---|
| elbow | 393.8 | 393.8 |
| nozzle | 390.5 | 390.5 |
