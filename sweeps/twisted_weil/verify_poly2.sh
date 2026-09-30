#!/bin/bash
# usage: verify_poly2.sh q g '[1,c1,...]' -> exhaustive search (even+odd models) for a genus-g hyperelliptic curve over F_q
S="$(cd "$(dirname "$0")" && pwd)"; REPO="$(cd "$S/../.." && pwd)"; mkdir -p "$S/logs"
[ -x "$S/hypsearch2" ] || cc -O2 -o "$S/hypsearch2" "$S/hypsearch2.c" || exit 1
q=$1; g=$2; P=$3
read N1 N2 <<< $(python3 -c "
c=$P; q=$q; e1=-c[1]; e2=c[2]; print(q+1-e1, q*q+1-(e1*e1-2*e2))")
tag=g${g}_$(echo "$q$P" | (md5 2>/dev/null || md5sum) | cut -c1-10)
$S/hypsearch2 $q $g $N1 $N2 > $S/logs/c2_$tag.txt 2> $S/logs/c2_$tag.err
MAGMA_STARTUP_FILE=/dev/null magma -b q:=$q IN:=$S/logs/c2_$tag.txt $S/candcheck2.m < /dev/null > $S/logs/c2_$tag.mout 2>&1
inT=$(grep -c -F -x "$P" $REPO/data/hypg${g}q$q.txt)
hit=$(grep -c -F "$P" $S/logs/c2_$tag.mout)
echo "VERIFY2 q=$q g=$g P=$P N1=$N1 N2=$N2 $(cat $S/logs/c2_$tag.err) intable=$inT found_by_search=$hit"
