#!/bin/bash
# rows for stage3_verdicts.txt from the 80-digit (and seeded) logs of the seventeen
cd /private/tmp/claude-501/-Users-assaferan-Documents-GitHub-ShimuraCurveALQuotients/2ccd0b0f-dcda-423c-b902-8aad70e4879f/scratchpad/nelold
fmt() { python3 -c "import sys; x=float(sys.argv[1]); print('%.2e'%x)" "$1" 2>/dev/null || echo "$1"; }
for it in 141_2:1,3,94,282 141_2:1,6,94,141 15_38:1,2,15,19,30,38,285,570 57_10:1,5,6,19,30,95,114,570 15_46:1,5,6,23,30,115,138,690 21_34:1,3,14,17,42,51,238,714 77_10:1,2,5,7,10,11,14,22,35,55,70,77,110,154,385,770 6_187:1,3,11,17,33,51,187,561 22_51:1,6,17,22,33,102,374,561 21_58:1,2,3,6,7,14,21,29,42,58,87,174,203,406,609,1218 39_34:1,2,3,6,13,17,26,34,39,51,78,102,221,442,663,1326 6_247:1,2,19,38,39,78,741,1482 6_253:1,2,23,33,46,66,759,1518 6_341:1,3,11,33,62,186,682,2046 6_455:1,6,7,10,13,15,42,70,78,91,105,130,195,546,910,1365 210_13:1,2,21,42,65,130,1365,2730 690_7:1,6,10,14,15,21,23,35,138,210,230,322,345,483,805,4830; do
  base=${it%%:*}; ws=${it#*:}; tag=$(echo $ws | tr ',' '-')
  for f in logs/s3o_${base}_W${tag}_c40_p80.log logs/s3o_${base}_W${tag}_c46_p80.log logs/s3o_${base}_W${tag}_c46_p80_seed1.log logs/s3o_${base}_W${tag}_c46_p120_seed1.log logs/s3o_${base}_W${tag}_c50_p120_seed1.log; do
    [ -f $f ] || continue
    grep -q 'STAGE3O DONE' $f || { echo "  $base W=$ws  $(basename $f): running"; continue; }
    cut=$(echo $f | sed 's/.*_c\([0-9]*\)_p.*/\1/'); pr=$(echo $f | sed 's/.*_p\([0-9]*\).*/\1/'); sd=""; case $f in *seed1*) sd=" seed 1";; esac
    g=$(grep -o 'genus [0-9]*' $f | head -1 | awk '{print $2}')
    v=$(grep VERDICT $f | sed 's/VERDICT: //'); s=$(fmt $(grep 'smallest singular' $f | awk '{print $4}'))
    gp=$(fmt $(grep 'largest ratio' $f | awk '{print $4}')); n=$(grep 'largest ratio' $f | awk '{print $6}')
    o=$(grep 'smallest singular' $f | sed 's/.*old columns //'); d=$(fmt $(grep 'largest deviation' $f | awk '{print $NF}'))
    t=$(grep 'lattice points total' $f | awk '{print int($6)}')
    printf "  %-7s g=%s  %-60s cut %s  %3s digits%s  %5ss  sv %s  gap %s after %s  W-check %s  old %s  %s\n" "$base" "$g" "{$ws}" "$cut" "$pr" "$sd" "$t" "$s" "$gp" "$n" "$d" "$o" "$v"
  done
done
