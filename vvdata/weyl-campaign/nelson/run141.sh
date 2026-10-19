#!/bin/bash
cd /private/tmp/claude-501/-Users-assaferan-Documents-GitHub-ShimuraCurveALQuotients/2ccd0b0f-dcda-423c-b902-8aad70e4879f/scratchpad/nelold
# wait for the retired worker's Magma (6_187 cut 46) to finish, record it, then run the seeded pair
until ! kill -0 62981 2>/dev/null; do sleep 20; done
f=logs/s3o_6_187_W1-3-11-17-33-51-187-561_c46_p80.log
echo "$(date '+%H:%M:%S') done  6_187 W=1,3,11,17,33,51,187,561 cut=46 prec=80 : $(grep VERDICT $f) $(grep 'smallest singular' $f) [orphan of retired worker]" >> logs/queue.log
for PR in 80 120; do
  log=logs/s3o_141_2_W1-6-94-141_c46_p${PR}_seed1.log
  echo "$(date '+%H:%M:%S') start 141_2 W=1,6,94,141 cut=46 prec=$PR SEED=1 npts=18 [coordinator's precision check]" >> logs/queue.log
  /Applications/Magma/magma -b DISC:=141 LEV:=2 WS:=1,6,94,141 PR:=$PR CUT:=46 NPTS:=18 CHECKW:=1 SEED:=1 stage3o_120.m < /dev/null > $log 2>&1
  echo "$(date '+%H:%M:%S') done  141_2 W=1,6,94,141 cut=46 prec=$PR SEED=1 : $(grep VERDICT $log) $(grep 'smallest singular' $log | cut -c1-60)" >> logs/queue.log
done
echo "$(date '+%H:%M:%S') seeded pair finished; restarting a worker" >> logs/queue.log
nohup ./runq2.sh > logs/w4.out 2>&1 &
