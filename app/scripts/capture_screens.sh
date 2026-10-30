#!/bin/sh
# Render every scenario in both themes: capture_screens.sh APP OUTDIR [scenario...]
# APP is the Kutta executable (macOS: Kutta.app/Contents/MacOS/Kutta).
set -eu
app=$1
out=$2
shift 2
mkdir -p "$out"
scenarios=${*:-$("$app" --list-scenarios)}
status=0
for s in $scenarios; do
  for theme in dark light; do
    if ! "$app" --scenario "$s" --theme "$theme" --size 1500x940 --screenshot "$out/$s-$theme.png"; then
      echo "FAILED: $s ($theme)" >&2
      status=1
    fi
  done
done
exit $status
