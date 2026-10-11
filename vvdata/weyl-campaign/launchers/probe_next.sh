#!/bin/bash
# 30-minute probes, one base per M, to read the first polytope rung each M-group asks for (polymake/nmzsolve.err).
cd ~/gymodels/composite/tree || exit 1
export NORMALIZ_BIN=/usr/bin/normaliz NMZ_TIMEOUT=120 BFCACHE=1
mkdir -p ../probesnext
for b in 6:65 6:77 6:85 6:91 6:95; do D=${b%%:*}; N=${b##*:}
  d=../probesnext/${D}_${N}; mkdir -p $d
  (cd ~/gymodels/composite/tree && setsid nohup timeout 30m magma -b D_s:=$D N_s:=$N VERB:=3 OUTDIR:=$d genmodels.m > $d/run.log 2> $d/run.err < /dev/null &)
done
sleep 2; ps -u $USER -o etime,args | grep "[m]agma.exe -b D_s" | grep -c probesnext; uptime | sed 's/.*average://'
