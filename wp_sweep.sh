#!/bin/bash
set -u
for b in "$@"; do
  D=${b%_*}; N=${b#*_}
  out="m0data/wp_${b}.txt"
  [ -s "$out" ] && grep -q DONE "$out" && { echo "skip $b"; continue; }
  /Applications/Magma/magma -b DD:=$D NN:=$N OUTNAME:="$out" wpoly_dump.m > /dev/null 2>&1 &
done
wait
echo "wpoly sweep done"
