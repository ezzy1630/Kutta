#!/bin/sh
# Render scenarios in both themes: capture_screens.sh APP OUTDIR [scenario...]
# APP is the Kutta executable (macOS: Kutta.app/Contents/MacOS/Kutta).
# Elbow scenarios use the real CLI and the Ansys mixing elbow when
# KUTTA_REAL_CLI and ELBOW_CASE are set (the .dat.h5 beside the case is the
# Fluent reference); car-* and bench always use the app's fixture CLI.
# param-run, cylinder-run and wing-* use the real CLI with PARAM_CASE,
# CYLINDER_CASE and WING_CASE when set, else the fixture's recordings.
set -eu
app=$1
out=$2
shift 2
mkdir -p "$out"
scenarios=${*:-$("$app" --list-scenarios)}
status=0
for s in $scenarios; do
  set --
  case $s in
    car-*|bench) ;;
    param-run|cylinder-run|wing-*)
      case $s in param-run) c=${PARAM_CASE:-} ;; cylinder-run) c=${CYLINDER_CASE:-} ;; *) c=${WING_CASE:-} ;; esac
      if [ -n "${KUTTA_REAL_CLI:-}" ] && [ -n "$c" ]; then set -- --cli "$KUTTA_REAL_CLI" "$c"; fi ;;
    *) if [ -n "${KUTTA_REAL_CLI:-}" ] && [ -n "${ELBOW_CASE:-}" ]; then set -- --cli "$KUTTA_REAL_CLI" "$ELBOW_CASE"; fi ;;
  esac
  for theme in dark light; do
    data=$(mktemp -d "${TMPDIR:-/tmp}/kutta-shot.XXXXXX")
    if ! "$app" --data-dir "$data" --scenario "$s" --theme "$theme" --size 1500x940 --screenshot "$out/$s-$theme.png" "$@"; then
      echo "FAILED: $s ($theme)" >&2
      status=1
    fi
    rm -rf "$data"
  done
done
exit $status
