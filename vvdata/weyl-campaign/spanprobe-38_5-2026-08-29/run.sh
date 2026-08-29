#!/bin/zsh
# run.sh <D> <N> <bump> <cap_seconds>
set -u
ROOT=/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients-spanprobe
OUT=/private/tmp/claude-501/-Users-assaferan-Documents-GitHub-ShimuraCurveALQuotients/c698c1d5-c594-43cf-9da5-ba723b4d6a52/scratchpad/spanprobe
export NORMALIZ_BIN=/Users/assaferan/Documents/GitHub/normaliz-3.11.1/normaliz
D=$1; N=$2; BUMP=$3; CAP=$4
LOG=$OUT/span_${D}_${N}_b${BUMP}.log
rm -f $LOG
t0=$(date +%s)
PROBE_BUMP=$BUMP magma -b DD:=$D NN:=$N $ROOT/spanprobe.m > $LOG 2>&1 &
pid=$!
while kill -0 $pid 2>/dev/null; do
  el=$(( $(date +%s) - t0 ))
  if [ $el -ge $CAP ]; then kill -9 $pid 2>/dev/null; wait $pid 2>/dev/null; echo "CAPPED ${el}s" >> $LOG; break; fi
  if grep -q "^DONE" $LOG 2>/dev/null; then sleep 2; kill -9 $pid 2>/dev/null; wait $pid 2>/dev/null; break; fi
  sleep 5
done
wait $pid 2>/dev/null
echo "ELAPSED $(( $(date +%s) - t0 ))s" >> $LOG
