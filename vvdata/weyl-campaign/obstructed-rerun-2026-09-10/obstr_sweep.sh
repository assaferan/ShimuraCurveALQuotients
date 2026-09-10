#!/bin/zsh
# Re-run every recorded Borcherds-obstructed base against CURRENT code.
# Verdicts in the sweep122 table were taken 2026-09-01/02, BEFORE the vx fix (d9b52d0, 09-05)
# touched the very stage that raises "Failed to find all Borcherds forms".
cd /Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign
export NORMALIZ_BIN=~/Documents/GitHub/normaliz-3.11.1/normaliz
OUT="$1"; LIST="$2"; CONC="${3:-5}"
run_one() {
  b=$1; D=${b%_*}; N=${b#*_}; log="$OUT/$b.log"; t0=$SECONDS
  PROBE_BUMP=0 magma -b DD:=$D NN:=$N vvdata/weyl-campaign/spanprobe.m < /dev/null > "$log" 2>&1
  el=$((SECONDS-t0))
  if grep -q "SPANPROBE RESULT .* SUCCESS" "$log"; then v="*** SUCCESS -- NO LONGER OBSTRUCTED ***"
  elif grep -q "Failed to find all Borcherds forms" "$log"; then v="obstructed (unchanged)"
  else v="OTHER: $(grep -iE 'error|assert' "$log" | head -1 | cut -c1-70)"; fi
  printf "%-9s %5ss  %s\n" "$b" "$el" "$v" >> "$OUT/VERDICTS.txt"
}
n=0
for b in $(cat "$LIST"); do
  run_one "$b" &
  n=$((n+1))
  if [ $((n % CONC)) -eq 0 ]; then wait; fi
done
wait
echo "SWEEP DONE" >> "$OUT/VERDICTS.txt"
