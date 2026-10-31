#!/bin/zsh
# Coupled before/after and SIMPLEC, on the V6 elbow (cold hybrid start) and the
# SST nozzle, CPU and Metal. Two criteria:
#  - fluent: `kutta solve` with the case's own Fluent residual criteria
#    (monitors/check-convergence switched on; the cases save it off), checked
#    every 20 iterations by the app's convergence monitor;
#  - strict: bench_convergence's V6 criterion, steady transport defects below
#    1e-6 and raw continuity below the saved Fluent run's final value, sampled
#    every 10 iterations.
# Baseline and candidate alternate within each repetition so that shared
# CPU/GPU load does not favour one. Run from the checkout root:
#   validation/results/coupled_convergence/measure.sh <baseline build> <candidate build> <outdir> [repetitions]
base=$1; cand=$2; out=$3; reps=${4:-2}; mkdir -p $out
fixtures=${KUTTA_FLUENT_EXAMPLES:-validation/fluent_examples}
typeset -A cases
cases=(elbow "$fixtures/mixing_elbow/mixing_elbow.cas.h5"
       nozzle "$fixtures/flprj_3891/nozzle_para_named_1.cffdb/nozzle_para_named_1-Solve/Base DP/nozzle_para_named_1.cas.h5")
print -r -- '{"v":1,"overrides":{"monitors/check-convergence":true}}' > $out/coupled.json
print -r -- '{"v":1,"overrides":{"monitors/check-convergence":true,"methods/pressure-velocity-coupling":"simplec"}}' > $out/simplec.json
for rep in {1..$reps}; do
  for backend in cpu metal; do
    for name in elbow nozzle; do
      # SIMPLEC is bitwise identical in both builds (parity.sh), so it runs once per repetition.
      order=(base cand simplec); (( rep % 2 )) || order=(simplec cand base)
      for side in $order; do
        bin=$cand; [[ $side == base ]] && bin=$base
        setup=$out/coupled.json; scheme=coupled; [[ $side == simplec ]] && { setup=$out/simplec.json; scheme=simplec; }
        tag=${name}_${backend}_${side}_r$rep
        echo "$(date +%T) $tag" >&2
        $bin/kutta solve "$cases[$name]" --approximate wfgc --backend $backend --iterations 1500 \
          --setup $setup --out $out/fluent_$tag --progress jsonl > $out/fluent_$tag.jsonl 2> $out/fluent_$tag.err
        $bin/bench_convergence --case $name --backend $backend --scheme $scheme --iterations 1500 \
          --interval 10 > $out/strict_$tag.log 2>&1
      done
    done
  done
done
