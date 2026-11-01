# Commands behind REPORT.md

Scripts as run (paths shortened). `main` binaries are a Release build of `e308968`; `cand` binaries are this branch.

## wing.sh

```sh
#!/bin/bash
# usage: wing.sh <bench> <mesh> <out> [extra args...]
bench=$1; mesh=$2; out=$3; shift 3
"$bench" --mesh "$mesh" --re 4.62e5 --downforce --moving-ground --tu 0.001 --ratio 1 --settle -1 --out "$out" "$@" > "$out.log" 2>&1
```

## ahmed.sh

```sh
#!/bin/bash
# usage: ahmed.sh <bench> <mesh> <out> [extra args...]
bench=$1; mesh=$2; out=$3; shift 3
"$bench" --mesh "$mesh" --re 2.78e6 --length 1.044 --area 0.112 --walls body --lift-z --tu 0.01 --window 1000 --settle -1 --out "$out" "$@" > "$out.log" 2>&1
```

## fsae2d.sh

```sh
#!/bin/bash
# usage: fsae2d.sh <tag> <bench> <outdir>
tag=$1; bench=$2; out=$3; mkdir -p $out
M=runs/fsae_aero/mesh
run() { name=$1; mesh=$2; shift 2; /usr/bin/time -l "$bench" --mesh $M/$mesh --out $out/$name --iterations 3000 --max-seconds 3600 --coupled "$@" > $out/$name.log 2>&1; }
run free zerihan_hfree_L2.msh --re 4.62e5 --downforce --tu 0.001 --ratio 1
run h0.448 zerihan_h0.448_L2.msh --re 4.62e5 --downforce --moving-ground --tu 0.001 --ratio 1
run h0.224 zerihan_h0.224_L2.msh --re 4.62e5 --downforce --moving-ground --tu 0.001 --ratio 1
run h0.134 zerihan_h0.134_L2.msh --re 4.62e5 --downforce --moving-ground --tu 0.001 --ratio 1
run h0.090 zerihan_h0.090_L2.msh --re 4.62e5 --downforce --moving-ground --tu 0.001 --ratio 1
run s4.76 s1223_a4.76_L2.msh --re 2e5 --tu 0.001 --ratio 1
run s10 s1223_a10.00_L2.msh --re 2e5 --tu 0.001 --ratio 1
echo finished > $out/DONE
```

## pair.sh

```sh
#!/bin/bash
cd <checkout>
D=runs/final; AM=runs/fsae_aero/mesh/ahmed35_L1.msh; WM=runs/fsae_aero/mesh/zerihan_h0.224_L1.msh
ah() { /usr/bin/time -l $1 --mesh $AM --re 2.78e6 --length 1.044 --area 0.112 --walls body --lift-z --tu 0.01 --window 1000 --settle -1 --coupled --iterations 2000 --out $D/ah_$2 > $D/ah_$2.log 2>&1; }
wi() { /usr/bin/time -l $1 --mesh $WM --re 4.62e5 --downforce --moving-ground --tu 0.001 --ratio 1 --settle -1 --coupled --iterations 3000 --out $D/w_$2 > $D/w_$2.log 2>&1; }
ah <main build>/bench_fsae_aero main & ah build/bench_fsae_aero cand &
wi <main build>/bench_fsae_aero main & wi build/bench_fsae_aero cand &
wait; echo done > $D/DONE
```

## internal.sh

```sh
#!/bin/zsh
cd <checkout>
F=$KUTTA_FLUENT_EXAMPLES; D=runs/internal; mkdir -p $D
nozzle="$F/flprj_3891/nozzle_para_named_1.cffdb/nozzle_para_named_1-Solve/Base DP/nozzle_para_named_1.cas.h5"
for side in main cand; do
  bin=<main build>/kutta; [[ $side == cand ]] && bin=<candidate build>/kutta
  /usr/bin/time -l $bin solve $F/mixing_elbow/mixing_elbow.cas.h5 --approximate wfgc --backend metal --iterations 3000 --out $D/elbow_$side > $D/elbow_$side.log 2>&1 &
  /usr/bin/time -l $bin solve "$nozzle" --approximate wfgc --backend metal --iterations 3000 --out $D/nozzle_$side > $D/nozzle_$side.log 2>&1 &
done
wait; echo done > $D/DONE
```

## slow.sh

```sh
#!/bin/bash
# usage: slow.sh <binary> <tag> <first> <last>
export KUTTA_FLUENT_EXAMPLES=$KUTTA_FLUENT_EXAMPLES
export KUTTA_FSAE_AERO_MESHES=runs/fsae_aero/mesh
out=scripts/slow/$2-$3
/usr/bin/time -l "$1" --test-suite=slow --first=$3 --last=$4 -r=xml --success > $out.xml 2> $out.err
echo "exit $?" >> $out.err
```

## Diagnosis runs

Built from `e308968` with `diag-hooks.patch`, run through `wing.sh` / `ahmed.sh`
with the environment shown (2000 iterations unless stated):

| run | environment |
| --- | --- |
| fturb | KUTTA_DIAG_FREEZE_TURB=600 |
| faphys | KUTTA_DIAG_FREEZE_APHYS=600 |
| cfl20, cfl5 | KUTTA_DIAG_CFL=20, 5 |
| lin3 | KUTTA_DIAG_LINTOL=1e-3 KUTTA_DIAG_LINMAX=100 (1000 iterations) |
| stdp | KUTTA_DIAG_STDP=1 |
| nolimit | KUTTA_DIAG_NOLIMIT=1 |
| fou | KUTTA_DIAG_FOU=1 |
| flim | KUTTA_DIAG_FREEZE_LIM=600 |
| lr01, lr002 | KUTTA_DIAG_FREEZE_LIM=600 KUTTA_DIAG_LIMRELAX=0.1, 0.02 (symmetric relaxation) |
| lr01_true | as lr01 with KUTTA_DIAG_TRUE_MONITOR=1 (monitor with the fresh limiter) |
| w_asym*, ah_asym* | KUTTA_DIAG_LIMASYM=1 KUTTA_DIAG_FREEZE_LIM=1 KUTTA_DIAG_LIMRELAX=r (falls at once) |
| ah_lr01, ah_s100, ah_s200 | symmetric relaxation 0.1 from iteration 1, 100, 200 |
| first_c | KUTTA_DIAG_NOFACEGRAD=1 KUTTA_DIAG_NOVELFACES=1 |
| nofacegrad, novf_c | KUTTA_DIAG_NOFACEGRAD=1; KUTTA_DIAG_NOVELFACES=1 |
| walknolim_c, ah_wnl | KUTTA_DIAG_WALKNOLIM=1 |
| switch_plain, switch_lr | SIMPLEC, KUTTA_DIAG_SWITCH=2000 (to Coupled 0.5/0.5), the second with symmetric relaxation 0.1 from 2001; 3500 iterations |
| sou_d | KUTTA_DIAG_DIR=… KUTTA_DIAG_DUMP_FROM=1000 KUTTA_DIAG_DUMP_TO=1011 (cell dumps; not archived, 37 MB) |

Reproduction runs (`repro/`, wing L1, `e308968`, `bench_fsae_aero --re 4.62e5
--downforce --moving-ground --tu 0.001 --ratio 1 --coupled --settle -1
--iterations 1500`, or 4000 for `w224L1_*`): `c_a025` adds `--urf-u 0.25 --urf-p
0.25`, `c_a09` adds `--urf-u 0.9 --urf-p 0.9`, `c_ts01` adds `--time-scale 0.1`,
and `c_check1` adds `--check 1`.

The held-limiter variants (`runs/win`) used the reverted commit `11189be` with
temporary `KUTTA_EXP_RATE`, `KUTTA_EXP_AGE` (snap window) and
`KUTTA_EXP_FRESH_AT` (fresh limiter from that iteration on) hooks.
