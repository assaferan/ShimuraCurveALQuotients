#!/bin/zsh
# queue.sh QUEUEFILE -- run the targets listed in QUEUEFILE one after another, each as a split
# evaluation (launch.sh), merged (merge.sh) and taken through fit / twist constants /
# certification (finish.sh, once per listed u).  Lines:  D N WS PR CUT NC NP NCH TAG [UMS]
# with UMS a comma-separated list of Atkin-Lehner indices for hypfit3's UM (empty = default).
# Progress goes to queue_<basename>.log; a target already carrying a merged ev_TAG.m is skipped.
#   nohup ./queue.sh queue1.txt < /dev/null > /dev/null 2>&1 &
set -u; setopt nonomatch
H=/private/tmp/claude-501/-Users-assaferan-Documents-GitHub-ShimuraCurveALQuotients/2ccd0b0f-dcda-423c-b902-8aad70e4879f/scratchpad/hyp16
Q=$1; L=$H/queue_$(basename $Q .txt).log
cd $H
echo "=== queue $Q started $(date)" >> $L
while read D N WS PR CUT NC NP NCH TAG UMS; do
  [ -z "$D" ] && continue
  [[ $D == \#* ]] && continue
  if [ -f $H/ev_$TAG.m ]; then echo "$(date) $TAG: ev_$TAG.m exists, skipping evaluation" >> $L
  else
    echo "$(date) $TAG: launching $NCH chunks (cut $CUT, $PR digits, $NC circles of $NP)" >> $L
    ./launch.sh $D $N $WS $PR $CUT $NC $NP $NCH $TAG >> $L 2>&1
    while true; do
      nd=$(grep -l 'HYPEVAL PART DONE' $H/ev_$TAG.part*.log 2>/dev/null | wc -l | tr -d ' ')
      [ "$nd" -eq "$NCH" ] && break
      if grep -q '^>>>\|Runtime error' $H/ev_$TAG.part*.log 2>/dev/null; then echo "$(date) $TAG: a chunk FAILED" >> $L; break; fi
      sleep 120
    done
    ./merge.sh $TAG >> $L 2>&1 || { echo "$(date) $TAG: merge failed, skipping" >> $L; continue; }
  fi
  if [ -z "${UMS:-}" ]; then
    echo "--- $(date) $TAG finish (default u)" >> $L; ./finish.sh $TAG >> $L 2>&1
  else
    for um in ${(s:,:)UMS}; do echo "--- $(date) $TAG finish u = w_$um" >> $L; ./finish.sh $TAG $um >> $L 2>&1; done
  fi
done < $Q
echo "=== queue $Q finished $(date)" >> $L
