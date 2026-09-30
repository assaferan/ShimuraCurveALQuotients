#!/bin/bash
# usage: verify_poly.sh q '[1,c1,...,c6]'  -> is there a genus-3 hyperelliptic curve over F_q with this Weil poly?
S="$(cd "$(dirname "$0")" && pwd)"; REPO="$(cd "$S/../.." && pwd)"; mkdir -p "$S/logs"
[ -x "$S/hypsearch" ] || cc -O2 -o "$S/hypsearch" "$S/hypsearch.c" || exit 1
q=$1; P=$2
read N1 N2 <<< $(python3 -c "
c=$P; q=$q; e1=-c[1]; e2=c[2]; print(q+1-e1, q*q+1-(e1*e1-2*e2))")
tag=$(echo "$q$P" | tr -d '[],-' | cut -c1-40)$(echo $P | (md5 2>/dev/null || md5sum) | cut -c1-6)
$S/hypsearch $q $N1 $N2 > $S/logs/cand_$tag.txt 2> $S/logs/cand_$tag.err
MAGMA_STARTUP_FILE=/dev/null magma -b q:=$q IN:=$S/logs/cand_$tag.txt $S/candcheck.m < /dev/null > $S/logs/cand_$tag.mout 2>&1
inT=$(grep -c -F -x "$P" $REPO/data/hypg3q$q.txt)
hit=$(grep -c -F "$P" $S/logs/cand_$tag.mout)
echo "VERIFY q=$q P=$P N1=$N1 N2=$N2 $(cat $S/logs/cand_$tag.err) intable=$inT found_by_search=$hit"
