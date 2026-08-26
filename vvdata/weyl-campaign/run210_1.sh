#!/bin/zsh
# 210_1 (composite-s kill-shot) end-to-end runbook.
# Prereq: $SOL145 = completed normaliz output for the (420, 145, 0) polytope.
# Usage: run210_1.sh <SOL145> <W0FILE> (both in nmzsolve/polymake_solution format)
set -e
SOL145=$1; W0=$2
CAMP=$(cd "$(dirname "$0")/../.." && pwd)   # campaign checkout root
CW=$CAMP/vvdata/weyl-campaign
SCI=${SCI:-$CW}                              # where to put logs
export NMZSOLVE=$CAMP/nmzsolve.py
export NORMALIZ_BIN=${NORMALIZ_BIN:-$HOME/Documents/GitHub/normaliz-3.11.1/normaliz}

# 1. seed the enumerated first rung
cp $SOL145 $CAMP/polymake/polymake_solution_420_145_0
# 2. install the weight-0 shift set; nmzsolve's t-shift fallback answers ALL
#    higher m=0 rung requests from the cache (thinned), whatever n the loop asks
cp $W0 $CAMP/polymake/tshift_w0_420.txt
# 3. mono stage (BorcherdsForms panel dump); Zero-side polytopes solve live via nmzsolve
cd $CAMP
magma -b DD:=210 NN:=1 $CW/monobase.m > $SCI/mono210_1.log 2>&1
head -4 $SCI/mono210_1.log
UK=$(grep -a "^USEDKEYS" $SCI/mono210_1.log | sed 's/USEDKEYS //; s/,$//; s/ //g')
echo "keys: $UK"
# 4. eta pool via the MITM enumerator (nd=24: use enum32m, NOT the old grid sweep)
python3 $CW/enum32m.py $SCI/mono210_1.log $SCI/epool_210_1.txt 8 7
# 5. eis32 with the cubic triple-loop DISABLED (nm at M=420 would make it lethal)
magma -b DD:=210 NN:=1 NOTRIPLES:=1 KEYS:=$UK EF:=$SCI/epool_210_1.txt $CW/eis32.m > $SCI/eis32e_210_1.out 2>&1
grep -a "SOLVE resid" $SCI/eis32e_210_1.out
echo "== continue with kernelrat/eisrat/genallgross as in chainbase.sh =="
