#!/bin/zsh
# merge.sh TAG  -- assemble ev_TAG.m from the header and the parts, in increasing JS order,
# refusing unless every part log says HYPEVAL PART DONE and the row count equals NC*(NP+1).
set -u
H=/private/tmp/claude-501/-Users-assaferan-Documents-GitHub-ShimuraCurveALQuotients/2ccd0b0f-dcda-423c-b902-8aad70e4879f/scratchpad/hyp16
TAG=$1
OUT=$H/ev_$TAG.m
read D N WS PR CUT NC NP NCH < $H/ev_$TAG.args
TOT=$(( NC * (NP + 1) ))
ndone=$(grep -l 'HYPEVAL PART DONE' $H/ev_$TAG.part*.log 2>/dev/null | wc -l | tr -d ' ')
if [ "$ndone" -ne "$NCH" ]; then echo "only $ndone of $NCH chunks finished"; exit 1; fi
parts=$(ls $OUT.p* | sed 's/.*\.p//' | sort -n | sed "s|^|$OUT.p|")
cat $OUT.hdr > $OUT
for p in ${(f)parts}; do cat $p >> $OUT; done
echo "];" >> $OUT
nrows=$(cat ${(f)parts} | wc -l | tr -d ' ')
if [ "$nrows" -ne "$TOT" ]; then echo "ROW COUNT $nrows != $TOT"; exit 1; fi
echo "merged $OUT: $nrows rows from $NCH parts"
grep -h 'evaluations:' $H/ev_$TAG.part*.log
