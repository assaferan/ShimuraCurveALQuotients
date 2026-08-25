#!/bin/zsh
# chainbase.sh <D> <N> -- full fourth-base pipeline at any base
set -e
D=$1; N=$2; TAG=${D}_${N}
SP=/private/tmp/claude-501/-Users-assaferan-Documents-GitHub-ShimuraCurveALQuotients/dbce22ad-34ec-4786-94b9-f6da57354d63/scratchpad; CAMP=/private/tmp/claude-501/-Users-assaferan-Documents-GitHub-ShimuraCurveALQuotients/e91ada1b-fe8a-4a7b-a418-93559b9c5c11/scratchpad/campaign; CW=/private/tmp/claude-501/-Users-assaferan-Documents-GitHub-ShimuraCurveALQuotients/e91ada1b-fe8a-4a7b-a418-93559b9c5c11/scratchpad/campaign/vvdata/weyl-campaign
cd $CAMP
magma -b DD:=$D NN:=$N $SP/sci/monobase.m > $SP/sci/mono$TAG.log 2>&1
head -4 $SP/sci/mono$TAG.log
UK=$(grep -a "^USEDKEYS" $SP/sci/mono$TAG.log | sed 's/USEDKEYS //; s/,$//; s/ //g')
echo "keys: $UK"
$SP/sci/venv/bin/python $SP/sci/enum32k.py $SP/sci/mono$TAG.log $SP/sci/epool_$TAG.txt 8 7
magma -b DD:=$D NN:=$N KEYS:=$UK EF:=$SP/sci/epool_$TAG.txt $CW/eis32.m > $SP/sci/eis32e_$TAG.out 2>&1
grep -a "SOLVE resid" $SP/sci/eis32e_$TAG.out
$SP/sci/venv/bin/python $CW/kernelrat.py $SP/sci/eis32e_$TAG.out > $SP/sci/kernelrat_$TAG.log 2>&1 || true
tail -2 $SP/sci/kernelrat_$TAG.log
$SP/sci/venv/bin/python $CW/eisrat.py $SP/sci/eis32e_$TAG.out > $SP/sci/eisrat_$TAG.log 2>&1 || true
tail -3 $SP/sci/eisrat_$TAG.log
$SP/sci/venv/bin/python $SP/sci/genallgross.py $SP/sci/mono$TAG.log $SP/sci/eisrat_$TAG.log $SP/sci/allgross$TAG.m
magma -b $SP/sci/allgross$TAG.m > $SP/sci/allgross_$TAG.log 2>&1
tail -30 $SP/sci/allgross_$TAG.log
