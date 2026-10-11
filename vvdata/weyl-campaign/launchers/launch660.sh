#!/bin/bash
# The nine M = 660 bases, once the (660, 231, 0) rung is installed. Usage: bash launch660.sh <tree> <logdir> <outdir> <max parallel>
# Expects polymake/polymake_solution_660_231_0 and polymake/tshift_w0_660.txt in <tree>; no git pull while jobs run.
tree=$1; logs=$2; out=$3; maxp=${4:-2}
cd "$tree" || exit 1
[ -f polymake/polymake_solution_660_231_0 ] || { echo "rung missing"; exit 1; }
[ -f polymake/tshift_w0_660.txt ] || { echo "shift set missing"; exit 1; }
mkdir -p "$logs" "$out"
export NORMALIZ_BIN=/usr/bin/normaliz BFPROGRESS=1 BFCACHE=1 NMZ_TIMEOUT=7200
for b in 6:55 10:33 15:11 15:22 22:15 33:5 33:10 55:3 330:1; do
  D=${b%%:*}; N=${b##*:}
  while [ "$(ps -u $USER -o command | grep -c '[m]agma.exe -b D_s')" -ge "$maxp" ]; do sleep 300; done
  echo "$(date '+%F %T') launching $D $N" >> "$logs/queue660.log"
  timeout 96h magma -b D_s:=$D N_s:=$N VERB:=3 OUTDIR:="$out" genmodels.m > "$logs/${D}_${N}.log" 2> "$logs/${D}_${N}.err" < /dev/null &
  sleep 60
done
wait
echo "$(date '+%F %T') queue done" >> "$logs/queue660.log"
