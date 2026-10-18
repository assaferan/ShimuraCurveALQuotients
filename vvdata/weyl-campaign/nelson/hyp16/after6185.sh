#!/bin/zsh
# When the eight extra circles (6_185_c46x) and queue A's standard pass (6_185_c46) have both
# finished, combine the two dumps and take the 22-circle set through fit / twist constants /
# certification.  Log: after6185.log.   nohup ./after6185.sh < /dev/null > /dev/null 2>&1 &
set -u; setopt nonomatch
H=/private/tmp/claude-501/-Users-assaferan-Documents-GitHub-ShimuraCurveALQuotients/2ccd0b0f-dcda-423c-b902-8aad70e4879f/scratchpad/hyp16
cd $H; L=$H/after6185.log
echo "started $(date)" >> $L
while true; do
  nd=$(grep -l 'HYPEVAL PART DONE' $H/ev_6_185_c46x.part*.log 2>/dev/null | wc -l | tr -d ' ')
  [ "$nd" -eq 2 ] && break
  grep -q '^>>>\|Runtime error' $H/ev_6_185_c46x.part*.log 2>/dev/null && { echo "extra-circle chunk FAILED $(date)" >> $L; exit 1; }
  sleep 120
done
./merge.sh 6_185_c46x >> $L 2>&1 || exit 1
echo "extra circles merged $(date)" >> $L
# queue A writes ev_6_185_c46.m by merge.sh and then runs finish.sh on it; wait for that fit log
until [ -f $H/ev_6_185_c46.m ] && [ -f $H/fit_6_185_c46.log ] && grep -q 'HYPFIT3 DONE\|Runtime error\|^>>>' $H/fit_6_185_c46.log; do sleep 120; done
sleep 30
/Applications/Magma/magma -b IN1:=$H/ev_6_185_c46.m IN2:=$H/ev_6_185_c46x.m OUT:=$H/ev_6_185_c46all.m $H/combine.m < /dev/null >> $L 2>&1
cp $H/ev_6_185_c46.args $H/ev_6_185_c46all.args
echo "--- combined fit $(date)" >> $L
./finish.sh 6_185_c46all >> $L 2>&1
echo "finished $(date)" >> $L
