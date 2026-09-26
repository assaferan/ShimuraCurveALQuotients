#!/usr/bin/env bash
# Twisted-trace sweep, two phases, each parallel (xargs -P) with GNU timeout per level:
#   1. pending levels (levels_pending.txt, curves_pending.txt), Weil mode q < 4g^2  -> out/D_N.out
#   2. supplement: report.py recomputes per-curve coverage and writes levels_supplement.txt
#      ("D N PMIN V3ONLY") + curves_supplement.txt.  V3ONLY=0 rows run all ops at primes p >= PMIN
#      -> out/D_N.supp.out (fills the p in [61, 4g^2) gap of the PB=59 local run for g >= 4).
#      V3ONLY=1 rows run only the V3 ops of curves with 9 in W at p = 2 mod 3, p >= PMIN
#      -> out/D_N.v3.out (pre-v3all runs used V3 ops at p = 1 mod 3 only).
#   JOBS=8 TIMEOUT=48h ./run.sh            (from anywhere; paths are resolved from this file)
# Optional: PHASE=1|2|all (default all), LEVELS=/CURVES= (override phase-1 lists), SUPPLEVELS=,
#           PB=<prime cap, default none>.
# Restartable: a level whose output file has a DONE line is skipped; a partial output (killed/timed
# out) is kept as *.partial-<time>.out (report.py still reads its finished RES lines) and rerun.
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
REL="sweeps/twisted_trace"
JOBS="${JOBS:-8}"; TIMEOUT="${TIMEOUT:-48h}"; PB="${PB:-0}"; PHASE="${PHASE:-all}"
LEVELS="${LEVELS:-$HERE/levels_pending.txt}"
CURVES="${CURVES:-$REL/curves_pending.txt}"
SUPPLEVELS="${SUPPLEVELS:-$HERE/levels_supplement.txt}"   # override only for testing a subset
mkdir -p "$HERE/out" "$HERE/logs"

if command -v timeout >/dev/null; then TO=timeout
elif command -v gtimeout >/dev/null; then TO=gtimeout
else echo "ERROR: GNU timeout not found (Linux coreutils; 'brew install coreutils' on macOS)" >&2; exit 1; fi
command -v magma >/dev/null || { echo "ERROR: magma not on PATH" >&2; exit 1; }
command -v python3 >/dev/null || { echo "ERROR: python3 not on PATH" >&2; exit 1; }

export REPO REL PB TIMEOUT TO CURVES
# one <curves-file> <suffix> D N [PMIN] [V3ONLY]
one() {
  C="$1"; SFX="$2"; D="$3"; N="$4"; PM="${5:-0}"; V3="${6:-0}"
  tag="${D}_${N}${SFX}"; o="$REPO/$REL/out/$tag.out"; lg="$REPO/$REL/logs/$tag.log"
  if [ -f "$o" ] && grep -q '^DONE ' "$o"; then echo "skip  $tag (done)"; return 0; fi
  [ -f "$o" ] && mv "$o" "$REPO/$REL/out/$tag.partial-$(date +%s).out"
  echo "start $tag $(date '+%F %T')"
  cd "$REPO" && "$TO" "$TIMEOUT" magma -b D:="$D" N:="$N" PB:="$PB" PMIN:="$PM" V3ONLY:="$V3" IN:="$C" \
      OUT:="$REL/out/$tag.out" "$REL/twist.m" < /dev/null > "$lg" 2>&1
  rc=$?
  if grep -q '^DONE ' "$o" 2>/dev/null; then echo "done  $tag $(date '+%F %T')"
  else echo "FAIL  $tag rc=$rc (124 = timeout; see $lg)"; fi
}
export -f one

ndone() {  # ndone <levels file> <suffix>
  local ok=0 tot=0
  while read -r D N _; do
    tot=$((tot+1)); grep -qs '^DONE ' "$HERE/out/${D}_${N}$2.out" && ok=$((ok+1))
  done < <(grep -E '^[0-9]+ [0-9]+' "$1")
  echo "$ok/$tot"
}

if [ "$PHASE" = all ] || [ "$PHASE" = 1 ]; then
  echo "== phase 1: pending levels (Weil range)"
  grep -E '^[0-9]+ [0-9]+' "$LEVELS" | awk '{print $1, $2}' |
    xargs -P "$JOBS" -n 2 bash -c 'one "$CURVES" "" "$0" "$1"'
fi
if [ "$PHASE" = all ] || [ "$PHASE" = 2 ]; then
  echo "== phase 2: supplement levels (coverage gaps, from report.py)"
  python3 "$HERE/report.py" > "$HERE/logs/report_before_supp.txt"
  grep -E '^[0-9]+ [0-9]+ [0-9]+ [01]' "$SUPPLEVELS" |
    xargs -P "$JOBS" -n 4 bash -c 'one "$REL/curves_supplement.txt" "$([ "$3" = 1 ] && echo .v3 || echo .supp)" "$0" "$1" "$2" "$3"'
fi

python3 "$HERE/report.py" > "$HERE/logs/report_final.txt"
inc=$(sed -n 's/^=== 6. INCOMPLETE COVERAGE: \([0-9]*\) .*/\1/p' "$HERE/logs/report_final.txt")
ro=$(sed -n 's/.*RULED OUT \([0-9]*\) (D>1.*/\1/p' "$HERE/logs/report_final.txt")
echo "STATUS: pending $(ndone "$LEVELS" "") levels DONE; supplement list now $(grep -c . "$HERE/levels_supplement.txt") levels; incomplete-coverage curves ${inc:-?}; U ruled out ${ro:-?}. Full report: $REL/logs/report_final.txt (rerun ./run.sh to resume)"
