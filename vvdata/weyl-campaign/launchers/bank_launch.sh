#!/bin/bash
# Enumerate the first m = 0 rung of a level, plus its weight-0 shift set for the t-shift fallback.
#   bank_launch.sh <tree> <M> <n>
tree=$1; M=$2; n=$3
d=~/bank$M; mkdir -p $d/logs
cd "$tree" || exit 1
export NORMALIZ_BIN=/usr/bin/normaliz
NMZ_TIMEOUT=864000 setsid nohup timeout 10d python3 nmzsolve.py $M $n 0 $d/polymake_solution_${M}_${n}_0 > $d/logs/nmz_${M}_${n}_0.log 2>&1 < /dev/null &
NMZ_TIMEOUT=86400 setsid nohup timeout 2d python3 nmzsolve.py $M $n 0 $d/w0_${M}_${n} 0 1 0 > $d/logs/nmz_${M}_w0.log 2>&1 < /dev/null &
sleep 2
ps -u $USER -o etime,args | grep "[n]mzsolve.py $M" | cut -c1-80
