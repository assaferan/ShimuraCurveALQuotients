#!/bin/zsh
# finish.sh TAG [UM]  -- after merge.sh TAG: hypfit3 -> hypfinish3 -> certgen, each stage read
# from the previous stage's log (the invariants as printed, kappa and mu as printed), stopping
# at the first stage that does not recognise or does not determine.  Logs: fit_TAG[_uUM].log,
# fin_TAG[_uUM].log, cert_TAG[_uUM].log.
set -u; setopt nonomatch
H=/private/tmp/claude-501/-Users-assaferan-Documents-GitHub-ShimuraCurveALQuotients/2ccd0b0f-dcda-423c-b902-8aad70e4879f/scratchpad/hyp16
TAG=$1; UMARG=${2:-}
read D N WS PR CUT NC NP NCH < $H/ev_$TAG.args
SUF=""; UMOPT=""
if [ -n "$UMARG" ]; then SUF="_u$UMARG"; UMOPT="UM:=$UMARG"; fi
cd $H
/Applications/Magma/magma -b IN:=$H/ev_$TAG.m $UMOPT $H/hypfit3.m < /dev/null > $H/fit_$TAG$SUF.log 2>&1
grep -q 'HYPFIT3 DONE' $H/fit_$TAG$SUF.log || { echo "hypfit3 did not finish; see fit_$TAG$SUF.log"; tail -5 $H/fit_$TAG$SUF.log; exit 1; }
UM=$(grep -o 'Atkin-Lehner involution w_[0-9]*' $H/fit_$TAG$SUF.log | head -1 | sed 's/.*w_//')
AA=$(grep -E '^   A = -?[0-9]+(/[0-9]+)?$' $H/fit_$TAG$SUF.log | sed 's/.*= //')
BB=$(grep -E '^   B = -?[0-9]+(/[0-9]+)?$' $H/fit_$TAG$SUF.log | sed 's/.*= //')
CC0=$(grep -E '^   C = -?[0-9]+(/[0-9]+)?$' $H/fit_$TAG$SUF.log | sed 's/.*= //')
echo "fit: u = w_$UM, A = ${AA:-NOT RECOGNISED}, B = ${BB:-NOT RECOGNISED}, C = ${CC0:-NOT RECOGNISED}"
grep -E 'WORST PREDICTION|conic residual|imaginary parts|recognising' $H/fit_$TAG$SUF.log
if [ -z "$AA" ] || [ -z "$BB" ] || [ -z "$CC0" ]; then echo "not all three invariants recognised; stopping"; exit 2; fi
/Applications/Magma/magma -b DISC:=$D LEV:=$N WS:=$WS UM:=$UM AA:=$AA BB:=$BB CC0:=$CC0 $H/hypfinish3.m < /dev/null > $H/fin_$TAG$SUF.log 2>&1
grep -q 'HYPFINISH3 DONE' $H/fin_$TAG$SUF.log || { echo "hypfinish3 did not finish; see fin_$TAG$SUF.log"; tail -5 $H/fin_$TAG$SUF.log; exit 3; }
KAP=$(grep -E '^kappa = ' $H/fin_$TAG$SUF.log | sed 's/kappa = //')
MU=$(grep -E '^mu = ' $H/fin_$TAG$SUF.log | sed 's/mu = //')
echo "finish: kappa = $KAP, mu = $MU"
grep -E 'kernel dimension|VERIFIED|MODEL|Hyperelliptic Curve' $H/fin_$TAG$SUF.log
/Applications/Magma/magma -b DISC:=$D LEV:=$N WS:=$WS UM:=$UM AA:=$AA BB:=$BB CC0:=$CC0 KAP:=$KAP MU:=$MU $H/certgen.m < /dev/null > $H/cert_$TAG$SUF.log 2>&1
grep -E 'CERTGEN|DISAGREE|AGREE|FAIL|IDENTICAL|DIFFER|splitting|mismatches|Riemann' $H/cert_$TAG$SUF.log
