#!/bin/zsh
# launch.sh D N WS PR CUT NC NP NCHUNK TAG
#   e.g.  ./launch.sh 55 7 1,35,55,77 40 40 12 16 4 55_7_c40
# Splits the NC*(NP+1) evaluations of hypeval.m into NCHUNK processes of hypeval_part.m.
# Output: ev_TAG.m.hdr, ev_TAG.m.p<JS>, logs ev_TAG.part<k>.log; merge with merge.sh TAG.
set -u; setopt nonomatch
H=/private/tmp/claude-501/-Users-assaferan-Documents-GitHub-ShimuraCurveALQuotients/2ccd0b0f-dcda-423c-b902-8aad70e4879f/scratchpad/hyp16
D=$1; N=$2; WS=$3; PR=$4; CUT=$5; NC=$6; NP=$7; NCH=$8; TAG=$9; CTRF=${10:-}
CTROPT=""; if [ -n "$CTRF" ]; then CTROPT="CTRF:=$H/$CTRF"; fi
TOT=$(( NC * (NP + 1) ))
OUT=$H/ev_$TAG.m
rm -f $OUT.hdr $OUT.p*
echo "$D $N $WS $PR $CUT $NC $NP $NCH" > $H/ev_$TAG.args
SETUP=$H/setup_${D}_$N.m
if [ ! -f $SETUP ]; then
  /Applications/Magma/magma -b DISC:=$D LEV:=$N OUT:=$SETUP $H/hypsetup.m < /dev/null > $H/setup_${D}_$N.log 2>&1
  grep -q 'HYPSETUP DONE' $H/setup_${D}_$N.log || { echo "hypsetup failed"; cat $H/setup_${D}_$N.log; exit 1; }
fi
echo "Eichler order: $(grep OBAS $SETUP)"
k=0; js=1
while [ $js -le $TOT ]; do
  k=$((k+1))
  je=$(( js + TOT / NCH - 1 ))
  if [ $k -eq $NCH ]; then je=$TOT; fi
  if [ $je -gt $TOT ]; then je=$TOT; fi
  cd $H && nohup /Applications/Magma/magma -b DISC:=$D LEV:=$N WS:=$WS PR:=$PR CUT:=$CUT NC:=$NC NP:=$NP RHO:=300 \
       JS:=$js JE:=$je OUT:=$OUT SETUP:=$SETUP $CTROPT $H/hypeval_part.m < /dev/null > $H/ev_$TAG.part$k.log 2>&1 &
  echo "chunk $k: points $js..$je  pid $!"
  js=$(( je + 1 ))
done
