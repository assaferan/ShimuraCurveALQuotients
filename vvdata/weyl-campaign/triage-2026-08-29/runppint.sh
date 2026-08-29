#!/bin/zsh
# Triage driver: run ppint.m on a list of bases, one at a time, with a wall-clock cap.
# Usage: runppint.sh <tag> <cap_seconds> D_N D_N ...
# Writes $OUT/ppint_D_N.log per base and appends a one-line verdict to $OUT/verdicts_<tag>.txt.
# macOS has no timeout(1), and `magma -b` HANGS after a runtime error rather than exiting,
# so each run is backgrounded and reaped explicitly.
set -u
# Run from the MAIN worktree, now that a924e1a merged origin/main into tier1-models: it has
# the Normaliz backend (PR #40) and NOT the DO-NOT-MERGE divisor prototype (95bd502), so it
# is the PRODUCTION pipeline.  Do NOT run this from the campaign worktree -- that branch
# carries the prototype, which rewrites the divisor search and (by design) turns 33_2
# integral, so its verdicts describe a pipeline that exists nowhere else.
# ppint.m itself lives only on the campaign branch, but uses only package intrinsics, so
# passing its path while standing in main resolves AttachSpec against the main package.
ROOT=/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients
SCRIPT=/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients-campaign/vvdata/weyl-campaign/ppint.m
export NORMALIZ_BIN=/Users/assaferan/Documents/GitHub/normaliz-3.11.1/normaliz
OUT=/private/tmp/claude-501/-Users-assaferan-Documents-GitHub-ShimuraCurveALQuotients/bb810bbc-25d6-42db-96a9-f66f444a30c7/scratchpad/ppint
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
    # magma hangs at 0% CPU after a runtime error: the log is complete, so reap on content
    if grep -q "^DONE" $LOG 2>/dev/null; then sleep 2; kill -9 $pid 2>/dev/null; wait $pid 2>/dev/null; break; fi
    sleep 5
  done
  wait $pid 2>/dev/null
  el=$(( $(date +%s) - t0 ))
  bv=$(grep "^BASEVERD" $LOG 2>/dev/null | tail -1)
  if [ -n "$bv" ]; then
    echo "${D}_${N} ${el}s $bv" >> $V
  elif grep -qi "error" $LOG 2>/dev/null; then
    echo "${D}_${N} ${el}s ERROR: $(grep -i 'error' $LOG | head -1 | cut -c1-140)" >> $V
  else
    echo "${D}_${N} ${el}s TIMEOUT_OR_NOOUTPUT" >> $V
  fi
done
echo "BATCH_${TAG}_COMPLETE" >> $V
