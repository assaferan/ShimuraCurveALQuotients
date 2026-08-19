#!/bin/bash
# Drive the odd-N m=0 oracle over a list of bases, one Magma process each (Magma has no usable
# internal time cap, and per-base isolation keeps one bad base from killing the sweep).
#
# usage: ./oracle_sweep.sh <outdir> <jobs> <D_N> [<D_N> ...]
#   e.g. ./oracle_sweep.sh out 2 6_5 6_7 10_3 10_7
#
# Each base writes out/<D>_<N>.log (full stdout, including the [RCE] blocks the fit consumes).
set -u
OUT=${1:?outdir}; JOBS=${2:?jobs}; shift 2
mkdir -p "$OUT"
MAGMA=/Applications/Magma/magma
TIMEOUT=3600

run_one() {
  local base=$1
  local D=${base%_*} N=${base#*_}
  local log="$OUT/${D}_${N}.log"
  if [ -s "$log" ] && grep -q "^\[ORACLE\] $D $N elapsed" "$log"; then
    echo "  skip $base (already done)"; return
  fi
  echo "  start $base"
  # `timeout` is not on macOS by default; emulate with a watchdog
  "$MAGMA" -b DD:=$D NN:=$N oracle_one.m > "$log" 2>&1 &
  local pid=$!
  ( sleep $TIMEOUT; kill -9 $pid 2>/dev/null ) & local wd=$!
  wait $pid 2>/dev/null
  kill -9 $wd 2>/dev/null
  if grep -q "^\[ORACLE\] $D $N elapsed" "$log"; then
    echo "  done  $base -- $(grep -c '^\[RCE\] cover' "$log") cover block(s)"
  else
    echo "  TIMEOUT/DIED $base"
  fi
}

i=0
for base in "$@"; do
  run_one "$base" &
  i=$((i+1))
  if [ $((i % JOBS)) -eq 0 ]; then wait; fi
done
wait
echo "sweep complete; logs in $OUT/"
