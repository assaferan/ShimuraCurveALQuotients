#!/bin/zsh
# Wave-4 triage driver: re-run ppint.m against the WeaklyHolomorphicBasis speedup.
#
#   retriage.sh <tag> <cap_seconds> D_N D_N ...
#
# Differs from runppint.sh (wave 3) in exactly one way that matters: ROOT points at the
# whbasis-speedup worktree.  Everything else -- the reap-on-DONE logic, the wall-clock cap, the
# one-line-per-base verdict file -- is unchanged, so wave-4 verdicts are directly comparable to
# wave-3's.
#
# WHY RE-RUN.  Wave 3 left 18 bases TIMED OUT at a 2400 s cap (unmeasured, NOT failing) and
# never started 122 more.  End-to-end at X0^38(5) went 856 s -> 30 s on this branch with an
# unchanged verdict, so most of those should now complete comfortably inside the cap.
#
# PRECONDITION: only run this once CI is green on whbasis-speedup.  The point of the exercise is
# to fill in unmeasured bases, which is worthless if the pipeline underneath is not trusted.
#
# NOTE macOS has no timeout(1), and `magma -b` HANGS after a runtime error rather than exiting,
# so each run is backgrounded and reaped explicitly on the DONE marker.
set -u
ROOT=/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients-whspeed
SCRIPT=/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients-campaign/vvdata/weyl-campaign/ppint.m
export NORMALIZ_BIN=/Users/assaferan/Documents/GitHub/normaliz-3.11.1/normaliz
OUT=${RETRIAGE_OUT:-$HOME/Documents/GitHub/ShimuraCurveALQuotients-campaign/vvdata/weyl-campaign/triage-wave4}
mkdir -p $OUT
TAG=$1; CAP=$2; shift 2
V=$OUT/verdicts_$TAG.txt
cd $ROOT
for b in "$@"; do
  D=${b%_*}; N=${b#*_}
  LOG=$OUT/ppint_${D}_${N}.log
  rm -f $LOG
  t0=$(date +%s)
  magma -b DD:=$D NN:=$N $SCRIPT > $LOG 2>&1 &
  pid=$!
  while kill -0 $pid 2>/dev/null; do
    now=$(date +%s); el=$((now-t0))
    if [ $el -ge $CAP ]; then kill -9 $pid 2>/dev/null; wait $pid 2>/dev/null; break; fi
    if grep -q "^DONE" $LOG 2>/dev/null; then sleep 2; kill -9 $pid 2>/dev/null; wait $pid 2>/dev/null; break; fi
    sleep 5
  done
  wait $pid 2>/dev/null
  el=$(( $(date +%s) - t0 ))
  bv=$(grep "^BASEVERD" $LOG 2>/dev/null | tail -1)
  if [ -n "$bv" ]; then
    echo "${D}_${N} ${el}s $bv" >> $V
  else
    err=$(grep -m1 -E "^(Runtime error|.*Runtime error)" $LOG 2>/dev/null | head -1)
    if [ -n "$err" ]; then
      echo "${D}_${N} ${el}s ERROR: $err" >> $V
    else
      echo "${D}_${N} ${el}s TIMEOUT_OR_NOOUTPUT" >> $V
    fi
  fi
done
echo "BATCH_${TAG}_COMPLETE" >> $V
