#!/bin/bash
# Run the bases of one level once its first rung is installed, at most <maxp> at a time.
#   launchM.sh <tree> <logdir> <outdir> <maxp> <M> <n> D:N [D:N ...]
# Expects polymake/polymake_solution_<M>_<n>_0 and polymake/tshift_w0_<M>.txt in <tree>; no git pull while jobs run.
tree=$1; logs=$2; out=$3; maxp=$4; M=$5; n=$6; shift 6
cd "$tree" || exit 1
[ -f polymake/polymake_solution_${M}_${n}_0 ] || { echo "rung missing"; exit 1; }
[ -f polymake/tshift_w0_${M}.txt ] || { echo "shift set missing"; exit 1; }
mkdir -p "$logs" "$out"
export NORMALIZ_BIN=/usr/bin/normaliz BFPROGRESS=1 BFCACHE=1 NMZ_TIMEOUT=7200
for b in "$@"; do
  D=${b%%:*}; N=${b##*:}
  while [ "$(ps -u $USER -o command | grep -c "[m]agma.exe -b D_s.*OUTDIR:=$out")" -ge "$maxp" ]; do sleep 300; done
  echo "$(date '+%F %T') launching $D $N" >> "$logs/queue.log"
  timeout 96h magma -b D_s:=$D N_s:=$N VERB:=3 OUTDIR:="$out" genmodels.m > "$logs/${D}_${N}.log" 2> "$logs/${D}_${N}.err" < /dev/null &
  sleep 60
done
wait
echo "$(date '+%F %T') queue done" >> "$logs/queue.log"
