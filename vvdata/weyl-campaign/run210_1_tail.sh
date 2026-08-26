#!/bin/zsh
# 210_1 chain tail: kernelrat -> eisrat -> genallgross -> allgross.
# Run AFTER run210_1.sh finishes eis32 (mono210_1.log + eis32e_210_1.out in CW).
set -e
SP=/private/tmp/claude-501/-Users-assaferan-Documents-GitHub-ShimuraCurveALQuotients/dbce22ad-34ec-4786-94b9-f6da57354d63/scratchpad
CAMP=/private/tmp/claude-501/-Users-assaferan-Documents-GitHub-ShimuraCurveALQuotients/e91ada1b-fe8a-4a7b-a418-93559b9c5c11/scratchpad/campaign
CW=$CAMP/vvdata/weyl-campaign
TAG=210_1
cd $CAMP
grep -a "SOLVE resid" $CW/eis32e_$TAG.out
$SP/sci/venv/bin/python $CW/kernelrat.py $CW/eis32e_$TAG.out > $CW/kernelrat_$TAG.log 2>&1 || true
tail -2 $CW/kernelrat_$TAG.log
$SP/sci/venv/bin/python $CW/eisrat.py $CW/eis32e_$TAG.out > $CW/eisrat_$TAG.log 2>&1 || true
tail -3 $CW/eisrat_$TAG.log
$SP/sci/venv/bin/python $SP/sci/genallgross.py $CW/mono${TAG}_flat.log $CW/eisrat_$TAG.log $CW/allgross$TAG.m
magma -b $CW/allgross$TAG.m > $CW/allgross_$TAG.log 2>&1
echo "=== 210_1 CHAIN TAIL DONE ==="
tail -30 $CW/allgross_$TAG.log
