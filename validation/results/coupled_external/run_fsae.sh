#!/bin/zsh
# Coupled vs SIMPLEC on the FSAE external-aero meshes (tools/fsae_aero), Metal.
# Both schemes stop on the bench's criterion (Fluent-scaled residuals <= 1e-3
# and forces settled) or their budget. Coupled's continuity scale is far smaller
# than SIMPLEC's, so summarize.py also finds the iteration where Coupled matched
# SIMPLEC's raw continuity (see REPORT.md). Two queues run in parallel; runs
# with a RESULT line already in their log are skipped.
#   validation/results/coupled_external/run_fsae.sh <bench_fsae_aero> <mesh dir> <out dir>
bench=$1; mesh=$2; out=$3; mkdir -p $out
wing=(--re 4.62e5 --downforce --tu 0.001 --ratio 1)
s1223=(--re 2e5 --tu 0.001 --ratio 1)
ahmed=(--re 2.78e6 --length 1.044 --area 0.112 --walls body --lift-z --tu 0.01 --ratio 1 --window 1000)
solve() { # tag mesh scheme budget args...
  local tag=$1 m=$2 scheme=$3 budget=$4; shift 4
  local flags=(); [[ $scheme == coupled ]] && flags=(--coupled)
  grep -q '^RESULT' $out/${tag}_$scheme.log 2>/dev/null && return  # finished earlier
  $bench --mesh $mesh/$m.msh --backend metal --check 50 --iterations $budget --max-seconds 5400 \
    --out $out/${tag}_$scheme $flags "$@" > $out/${tag}_$scheme.log 2>&1
}
queue_a() {
  for h in free 0.448 0.224; do
    ground=(); [[ $h != free ]] && ground=(--moving-ground)
    solve zerihan_$h zerihan_h${h}_L2 coupled 2000 $wing $ground
    solve zerihan_$h zerihan_h${h}_L2 simplec 6000 $wing $ground
  done
  solve ahmed35_L1 ahmed35_L1 coupled 3000 $ahmed
}
queue_b() {
  for h in 0.134 0.090; do
    solve zerihan_$h zerihan_h${h}_L2 coupled 2000 $wing --moving-ground
    solve zerihan_$h zerihan_h${h}_L2 simplec 6000 $wing --moving-ground
  done
  for a in 4.76 10.00; do
    solve s1223_$a s1223_a${a}_L2 coupled 2000 $s1223
    solve s1223_$a s1223_a${a}_L2 simplec 6000 $s1223
  done
  solve ahmed35_L1 ahmed35_L1 simplec 5000 $ahmed
}
queue_a & queue_b & wait
