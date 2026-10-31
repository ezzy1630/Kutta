#!/bin/zsh
# Well-conditioned Fluent cases, baseline against candidate: `kutta solve
# --iterations 40` with each binary on CPU and Metal. Each case runs its saved
# scheme (V6 elbow and nozzle: Coupled; elbow_param: SIMPLE) and the elbow
# also runs SIMPLEC. Compares fluxes.csv, forces.csv, loads.json, result.vtu
# and residuals.csv without its last (seconds_per_iteration) field.
#   validation/results/external_3d_robustness/parity.sh <baseline kutta> <candidate kutta> <outdir>
base=$1; cand=$2; out=$3; mkdir -p $out
fixtures=${KUTTA_FLUENT_EXAMPLES:-validation/fluent_examples}
print -r -- '{"v":1,"overrides":{"methods/pressure-velocity-coupling":"simplec"}}' > $out/simplec.json
nozzle="$fixtures/flprj_3891/nozzle_para_named_1.cffdb/nozzle_para_named_1-Solve/Base DP/nozzle_para_named_1.cas.h5"
typeset -A args
args=(elbow_coupled "$fixtures/mixing_elbow/mixing_elbow.cas.h5|--approximate|wfgc"
      elbow_simplec "$fixtures/mixing_elbow/mixing_elbow.cas.h5|--approximate|wfgc|--setup|$out/simplec.json"
      nozzle_coupled "$nozzle|--approximate|wfgc"
      elbow_param_simple "$fixtures/mixing_elbow/elbow_param.cas.h5|--approximate|compressibility,ke-options")
result=0
for backend in cpu metal; do
  for name in elbow_coupled elbow_simplec nozzle_coupled elbow_param_simple; do
    for side in base cand; do
      bin=$base; [[ $side == cand ]] && bin=$cand
      dir=$out/${name}_${backend}_$side
      $bin solve "${(@s:|:)args[$name]}" --backend $backend --iterations 40 --out $dir > $dir.log 2>&1 || { echo "FAILED $dir"; result=1; }
    done
    a=$out/${name}_${backend}_base; b=$out/${name}_${backend}_cand; same=1
    for f in fluxes.csv forces.csv loads.json result.vtu; do cmp -s $a/$f $b/$f || { same=0; echo "DIFF $name $backend $f"; }; done
    cmp -s <(sed -E 's/,[^,]*$//' $a/residuals.csv) <(sed -E 's/,[^,]*$//' $b/residuals.csv) \
      || { same=0; echo "DIFF $name $backend residuals.csv"; }
    echo "$name $backend $([[ $same == 1 ]] && echo identical || echo DIFFERENT)" \
      "$(grep -o 'Pressure-velocity coupling: [A-Za-z]*' $b.log)"
    (( same )) || result=1
  done
done
exit $result
