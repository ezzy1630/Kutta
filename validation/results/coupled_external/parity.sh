#!/bin/zsh
# SIMPLEC and SIMPLE must be bitwise unchanged by the Coupled Jacobian change.
# Runs `kutta solve --iterations 40` with both binaries on CPU and Metal:
#   - the V6 mixing elbow with SIMPLEC forced by --setup (--approximate wfgc);
#   - elbow_param, which runs SIMPLE (--approximate compressibility,ke-options).
# Compares fluxes.csv, forces.csv, loads.json, result.vtu, and residuals.csv
# without its last (seconds_per_iteration) column. Run from the checkout root:
#   validation/results/coupled_external/parity.sh <baseline kutta> <candidate kutta> <outdir>
base=$1; cand=$2; out=$3; mkdir -p $out
fixtures=${KUTTA_FLUENT_EXAMPLES:-validation/fluent_examples}
print -r -- '{"v":1,"overrides":{"methods/pressure-velocity-coupling":"simplec"}}' > $out/simplec.json
result=0
run() { # name backend side bin args...
  local name=$1 backend=$2 side=$3 bin=$4; shift 4
  local dir=$out/${name}_${backend}_$side
  $bin solve "$@" --backend $backend --iterations 40 --out $dir > $dir.log 2>&1 || { echo "FAILED $dir"; result=1; }
}
for backend in cpu metal; do
  for side in base cand; do
    bin=$base; [[ $side == cand ]] && bin=$cand
    run elbow_simplec $backend $side $bin "$fixtures/mixing_elbow/mixing_elbow.cas.h5" --approximate wfgc --setup $out/simplec.json
    run elbow_param_simple $backend $side $bin "$fixtures/mixing_elbow/elbow_param.cas.h5" --approximate compressibility,ke-options
  done
  for name in elbow_simplec elbow_param_simple; do
    a=$out/${name}_${backend}_base; b=$out/${name}_${backend}_cand
    same=1
    for f in fluxes.csv forces.csv loads.json result.vtu; do cmp -s $a/$f $b/$f || { same=0; echo "DIFF $name $backend $f"; }; done
    # The file opens with comment lines; the last field of each line is dropped.
    cmp -s <(sed -E 's/,[^,]*$//' $a/residuals.csv) <(sed -E 's/,[^,]*$//' $b/residuals.csv) \
      || { same=0; echo "DIFF $name $backend residuals.csv"; }
    echo "$name $backend $([[ $same == 1 ]] && echo identical || echo DIFFERENT)" \
      "$(grep -o 'Pressure-velocity coupling: [A-Za-z]*' $b.log)" \
      "sha256(result.vtu)=$(shasum -a 256 < $b/result.vtu | cut -c1-16)"
    (( same )) || result=1
  done
done
exit $result
