#!/usr/bin/env bash
# Run the twisted-trace test on every level in levels_pending.txt, in parallel.
#   JOBS=8 TIMEOUT=48h ./run.sh            (from anywhere; paths are resolved from this file)
# Optional: LEVELS=<file of "D N" lines>  PB=<prime bound, default 59>
#           CURVES=<curve list, repo-relative, default sweeps/twisted_trace/curves_pending.txt>
# Each level writes out/D_N.out (results) and logs/D_N.log (Magma stdout/stderr).
# Levels whose out/D_N.out already has a DONE line are skipped, so this is safe to restart.
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
REL="sweeps/twisted_trace"
JOBS="${JOBS:-8}"; TIMEOUT="${TIMEOUT:-48h}"; PB="${PB:-59}"
LEVELS="${LEVELS:-$HERE/levels_pending.txt}"
CURVES="${CURVES:-$REL/curves_pending.txt}"
mkdir -p "$HERE/out" "$HERE/logs"

if command -v timeout >/dev/null; then TO=timeout
elif command -v gtimeout >/dev/null; then TO=gtimeout
else echo "ERROR: GNU timeout not found (Linux coreutils; 'brew install coreutils' on macOS)" >&2; exit 1; fi
command -v magma >/dev/null || { echo "ERROR: magma not on PATH" >&2; exit 1; }

export REPO REL PB TIMEOUT TO CURVES
one() {
  D="$1"; N="$2"; o="$REPO/$REL/out/${D}_${N}.out"; lg="$REPO/$REL/logs/${D}_${N}.log"
  if [ -f "$o" ] && grep -q '^DONE ' "$o"; then echo "skip $D $N (done)"; return 0; fi
  echo "start $D $N $(date '+%F %T')"
  cd "$REPO" && "$TO" "$TIMEOUT" magma -b D:="$D" N:="$N" PB:="$PB" IN:="$CURVES" \
      OUT:="$REL/out/${D}_${N}.out" "$REL/twist.m" < /dev/null > "$lg" 2>&1
  rc=$?
  if grep -q '^DONE ' "$o" 2>/dev/null; then echo "done  $D $N $(date '+%F %T')"
  else echo "FAIL  $D $N rc=$rc (124 = timeout; see $lg)"; fi
}
export -f one

grep -E '^[0-9]+ [0-9]+' "$LEVELS" | xargs -P "$JOBS" -n 2 bash -c 'one "$0" "$1"'

tot=$(grep -cE '^[0-9]+ [0-9]+' "$LEVELS"); ok=0
while read -r D N _; do
  [ -f "$HERE/out/${D}_${N}.out" ] && grep -q '^DONE ' "$HERE/out/${D}_${N}.out" && ok=$((ok+1))
done < <(grep -E '^[0-9]+ [0-9]+' "$LEVELS")
echo "STATUS: $ok/$tot levels DONE, $((tot-ok)) not done (rerun ./run.sh to resume); then: python3 $REL/report.py"
