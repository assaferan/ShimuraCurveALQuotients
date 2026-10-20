#!/bin/bash
# Twisted Weil-polynomial sweep: one Magma job per (D,N) level, JOBS in parallel, from anywhere.
#   ./run.sh [levels file, default levels.txt] [JOBS, default 5]
# Each level writes out/D_N.out (log in logs/D_N.log); a level whose output ends in DONE is skipped.
# Curves come from IN (default sweeps/twisted_weil/all_in.txt, relative to the repo root).
# Runs without .magmarc (MAGMA_STARTUP_FILE=/dev/null), like the local run.
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
REL="sweeps/twisted_weil"
LEVELS="${1:-levels.txt}"; case "$LEVELS" in /*) ;; *) LEVELS="$HERE/$LEVELS";; esac
IN="${IN:-$REL/all_in.txt}"
mkdir -p "$HERE/out" "$HERE/logs"
export REPO REL IN
one() { D=$1; N=$2; o="$REPO/$REL/out/${D}_${N}.out"
  if [ -f "$o" ] && tail -1 "$o" | grep -q '^DONE'; then return; fi
  cd "$REPO" && MAGMA_STARTUP_FILE=/dev/null magma -b D:=$D N:=$N IN:="$IN" OUT:="$REL/out/${D}_${N}.out" \
    "$REL/tw.m" < /dev/null > "$REPO/$REL/logs/${D}_${N}.log" 2>&1; }
export -f one
grep -E '^[0-9]+ [0-9]+' "$LEVELS" | awk '{print $1, $2}' | xargs -P "${2:-5}" -L 1 bash -c 'one $0 $1'
