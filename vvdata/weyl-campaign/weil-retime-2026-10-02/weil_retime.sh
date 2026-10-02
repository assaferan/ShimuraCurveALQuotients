#!/bin/bash
# Weil-stage re-timing on lovelace: WeilPolynomial(X, p) per prime, main at 95b19e6, class-number
# tables mounted.  33 curves spanning genus 3-7, #W 1-64 and level 30-30030, including the five
# already timed (1071, 7296, 13029, 9755, 785).  Run from ~/gymodels/tree_main; logs in
# ~/gymodels/weil/logs/weil_main2_<id>.log (per-prime lines) and run_<id>.out (stdout).
cd ~/gymodels/tree_main || exit 1
export CLASS_GROUPS_FAST_DIR=/scratch/class-groups-fast
LOGS=$HOME/gymodels/weil/logs
mkdir -p $LOGS
IDS="94 2176 976 12610 1387 3853 6260 18378 136 2174 785 2122 979 2586 9755 3854 7060 13029 153 2108 582 977 1416 1434 9325 232 2325 1071 2273 1219 1436 346 7296 2420 1217 239"
echo "$IDS" | tr ' ' '\n' | xargs -P 8 -I{} sh -c 'timeout 12h magma -b id:={} tag:=main2 scratch:='"$LOGS"' weil_timing.m < /dev/null > '"$LOGS"'/run_{}.out 2>&1'
echo "RETIME_ALL_DONE" > $LOGS/ALL_DONE
