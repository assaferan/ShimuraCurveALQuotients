#!/bin/bash
# usage: runq.sh <prec> <cut> <checkw> item item ...   with item = D_N:w1,w2,...[:npts]
PR=$1; CUT=$2; CHK=$3; shift 3
cd /private/tmp/claude-501/-Users-assaferan-Documents-GitHub-ShimuraCurveALQuotients/2ccd0b0f-dcda-423c-b902-8aad70e4879f/scratchpad/nelold
for it in "$@"; do
  base=${it%%:*}; rest=${it#*:}; ws=${rest%%:*}; np=18
  case "$rest" in *:*) np=${rest##*:};; esac
  D=${base%_*}; N=${base#*_}
  tag=$(echo $ws | tr ',' '-')
  log=logs/s3o_${base}_W${tag}_c${CUT}_p${PR}.log
  echo "$(date '+%H:%M:%S') start $base W=$ws cut=$CUT prec=$PR npts=$np" >> logs/queue.log
  /Applications/Magma/magma -b DISC:=$D LEV:=$N WS:=$ws PR:=$PR CUT:=$CUT NPTS:=$np CHECKW:=$CHK stage3o.m < /dev/null > $log 2>&1
  echo "$(date '+%H:%M:%S') done  $base W=$ws cut=$CUT prec=$PR : $(grep VERDICT $log) $(grep 'smallest singular' $log)" >> logs/queue.log
done
