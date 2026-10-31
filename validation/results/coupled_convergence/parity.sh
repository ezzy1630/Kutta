#!/bin/zsh
# SIMPLEC and SIMPLE must be bitwise unchanged by the Coupled time step change.
# Runs `kutta solve --iterations 40` with each binary on the V6 elbow and the
# nozzle (SIMPLEC and SIMPLE forced by --setup), CPU and Metal, and compares
# fluxes.csv, forces.csv, loads.json and residuals.csv without its
# seconds_per_iteration column. Run from the checkout root:
#   validation/results/coupled_convergence/parity.sh <baseline kutta> <candidate kutta> <outdir>
base=$1; cand=$2; out=$3; mkdir -p $out
fixtures=${KUTTA_FLUENT_EXAMPLES:-validation/fluent_examples}
typeset -A cases
cases=(elbow "$fixtures/mixing_elbow/mixing_elbow.cas.h5"
       nozzle "$fixtures/flprj_3891/nozzle_para_named_1.cffdb/nozzle_para_named_1-Solve/Base DP/nozzle_para_named_1.cas.h5")
result=0
for scheme in simplec simple; do
  print -r -- "{\"v\":1,\"overrides\":{\"methods/pressure-velocity-coupling\":\"$scheme\"}}" > $out/$scheme.json
  for name in elbow nozzle; do
    for backend in cpu metal; do
      for side in base cand; do
        bin=$base; [[ $side == cand ]] && bin=$cand
        dir=$out/${name}_${scheme}_${backend}_$side
        $bin solve "$cases[$name]" --approximate wfgc --backend $backend --iterations 40 \
          --setup $out/$scheme.json --out $dir > $dir.log 2>&1 || { echo "FAILED $dir"; result=1; }
      done
      a=$out/${name}_${scheme}_${backend}_base; b=$out/${name}_${scheme}_${backend}_cand
      same=1
      for f in fluxes.csv forces.csv loads.json; do cmp -s $a/$f $b/$f || { same=0; echo "DIFF $f"; }; done
      cmp -s <(cut -d, -f1-$(( $(head -2 $a/residuals.csv | tail -1 | tr ',' '\n' | wc -l) - 1 )) $a/residuals.csv) \
             <(cut -d, -f1-$(( $(head -2 $b/residuals.csv | tail -1 | tr ',' '\n' | wc -l) - 1 )) $b/residuals.csv) \
        || { same=0; echo "DIFF residuals.csv"; }
      echo "$name $scheme $backend $([[ $same == 1 ]] && echo identical || echo DIFFERENT) sha256(fluxes)=$(shasum -a 256 < $b/fluxes.csv | cut -c1-16)"
      (( same )) || result=1
    done
  done
done
exit $result
