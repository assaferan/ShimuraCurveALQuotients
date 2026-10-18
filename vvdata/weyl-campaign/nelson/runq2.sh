#!/bin/bash
# claim-based worker: items17.txt lines "PR CUT CHK D_N:ws[:npts]"; a worker claims an item by mkdir and runs it.
cd /private/tmp/claude-501/-Users-assaferan-Documents-GitHub-ShimuraCurveALQuotients/2ccd0b0f-dcda-423c-b902-8aad70e4879f/scratchpad/nelold
mkdir -p claims logs
while read PR CUT CHK it; do
  [ -z "$it" ] && continue
  base=${it%%:*}; rest=${it#*:}; ws=${rest%%:*}; np=18
  case "$rest" in *:*) np=${rest##*:};; esac
  D=${base%_*}; N=${base#*_}; tag=$(echo $ws | tr ',' '-')
  key=${base}_W${tag}_c${CUT}_p${PR}
  mkdir claims/$key 2>/dev/null || continue
  log=logs/s3o_$key.log
  echo "$(date '+%H:%M:%S') start $base W=$ws cut=$CUT prec=$PR npts=$np [worker $$]" >> logs/queue.log
  /Applications/Magma/magma -b DISC:=$D LEV:=$N WS:=$ws PR:=$PR CUT:=$CUT NPTS:=$np CHECKW:=$CHK stage3o.m < /dev/null > $log 2>&1
  echo "$(date '+%H:%M:%S') done  $base W=$ws cut=$CUT prec=$PR : $(grep VERDICT $log) $(grep 'smallest singular' $log)" >> logs/queue.log
done < items17.txt
echo "$(date '+%H:%M:%S') worker $$ found nothing left" >> logs/queue.log
