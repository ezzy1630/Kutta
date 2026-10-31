#!/bin/zsh
# SIMPLEC vs Coupled on the bench cases, alternating which scheme runs first
# so shared-GPU drift does not favour one. Run from the checkout root with
# KUTTA_FLUENT_EXAMPLES set (the elbow) and the TMR grids prepared:
#   validation/results/coupled/run_comparison.sh build/bench_convergence <outdir>
# Coupled uses Fluent's pseudo-transient explicit relaxation 0.5/0.5; SIMPLEC
# uses each case's bench defaults. Peak RSS comes from /usr/bin/time -l.
bench=$1; out=$2; mkdir -p $out
i=0
run() { # name backend args...
  local name=$1 backend=$2; shift 2
  local order=(simplec coupled); (( i++ % 2 )) && order=(coupled simplec)
  for sc in $order; do
    local extra=()
    [[ $sc == coupled ]] && extra=(--momentum-relaxation .5 --pressure-relaxation .5)
    echo "$(date +%T) $name $backend $sc" >&2
    /usr/bin/time -l $bench --backend $backend --scheme $sc --max-seconds 900 "$@" $extra > $out/${name}_${backend}_${sc}.log 2>&1
    echo "EXIT $?" >> $out/${name}_${backend}_${sc}.log
  done
}
for backend in metal cpu; do
  run elbow $backend --case elbow --iterations 3000
  run plate-sst $backend --case plate --level 2 --iterations 3000
  run bump-sst $backend --case bump --level 2 --iterations 3000
  run plate-sa $backend --case plate --model sa --level 2 --iterations 3000
  run bump-sa $backend --case bump --model sa --level 2 --iterations 3000
  run channel $backend --case channel --iterations 1000
done
run channel-983k metal --case channel --nx 128 --ny 32 --nz 40 --iterations 200
