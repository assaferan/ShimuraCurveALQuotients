#!/bin/zsh
# 330_1 chain tail: kernelrat -> eisrat -> genallgross -> allgross.
# Run AFTER eis32k finishes (eis32k_330_1.out in CW).  The mono dump is the
# SYNTHETIC one (genallgross only reads its BASE line).
set -e
SP=/private/tmp/claude-501/-Users-assaferan-Documents-GitHub-ShimuraCurveALQuotients/dbce22ad-34ec-4786-94b9-f6da57354d63/scratchpad
CAMP=/private/tmp/claude-501/-Users-assaferan-Documents-GitHub-ShimuraCurveALQuotients/e91ada1b-fe8a-4a7b-a418-93559b9c5c11/scratchpad/campaign
CW=$CAMP/vvdata/weyl-campaign
TAG=330_1
cd $CAMP
grep -a "SOLVE resid" $CW/eis32k_$TAG.out
$SP/sci/venv/bin/python $CW/kernelrat.py $CW/eis32k_$TAG.out > $CW/kernelrat_$TAG.log 2>&1 || true
grep -a "rank" $CW/kernelrat_$TAG.log | tail -2
$SP/sci/venv/bin/python $CW/eisrat.py $CW/eis32k_$TAG.out > $CW/eisrat_$TAG.log 2>&1 || true
tail -3 $CW/eisrat_$TAG.log
$SP/sci/venv/bin/python $SP/sci/genallgross.py $CW/mono330_1_synth.log $CW/eisrat_$TAG.log $CW/allgross$TAG.m
magma -b $CW/allgross$TAG.m > $CW/allgross_$TAG.log 2>&1
echo "--- 330_1 CHAIN TAIL DONE ---"
tail -40 $CW/allgross_$TAG.log
