#!/bin/sh
# Collect the stage-3 runs into one table.  It prints what was measured and does NOT decide
# anything: a verdict needs two cutoffs and a human looking at the spectrum, for the reasons in
# stage3_verdicts.txt.
#
#   sh collect_stage3.sh <log-directory>
#
# Columns: the run, whether it finished, the verdict line the script printed, the smallest
# singular value (the residual), the SECOND singular value, and the gap with the number of
# relations.
#
# ⚠ NEXT-SV is literally the second singular value, nothing cleverer.  At genus 3 that is the
# first one above the gap and comparing it with the residual is the right thing to do.  AT GENUS
# 4 IT IS NOT: there are three relations, so the second singular value is itself one of the three
# near-zeros, and residual-against-NEXT-SV is a comparison of two zeros.  At genus 4 read the GAP
# column, which is taken after the correct number of values.
#
# How to read it, in order:
#   1. the control n_ctl_91w7_* must say NOT hyperelliptic with a residual near 1.8e-2.  If it
#      does not, nothing else in the batch is worth reading.
#   2. genus 3: one relation is hyperelliptic, none is not.  GENUS 4 (6_403, 2418_1, 3390_1):
#      THREE relations is hyperelliptic and one is not.
#   3. a hyperelliptic verdict needs the residual to FALL from cutoff 40 to 46 while the next
#      singular value stays put; a non-hyperelliptic one needs the residual to stay near 1e-2 and
#      not fall.  One cutoff is not a verdict.
#   4. a residual between 1e-2 and 1e-12 is to be investigated, not called.  Nothing measured so
#      far has landed there: real ones sit near 1e-2, numerical zeros below 1e-12.
#   5. if a residual stops responding to the cutoff, raise PR.  A real singular value does not
#      move with the precision (29 digits stable at the control); a numerical zero does (a factor
#      of 112 at X_0^1326(1)).
#
# A log is finished only when it ends STAGE3N DONE.  Magma buffers to file, so an empty log is a
# running job, not a stalled one, and a job killed for memory can leave 0 bytes looking exactly
# like one that never started.
set -u
DIR="${1:-.}"
printf '%-22s %-9s %-18s %-14s %-14s %s\n' RUN STATE VERDICT RESIDUAL NEXT-SV GAP
for f in "$DIR"/p_*.log "$DIR"/o_*.log "$DIR"/n_ctl_*.log "$DIR"/hi_*.log "$DIR"/nx_*.log; do
    [ -e "$f" ] || continue
    b=$(basename "$f" .log)
    if grep -q 'STAGE3N\? DONE' "$f" 2>/dev/null; then
        state=done
        v=$(grep '^VERDICT' "$f" | sed 's/VERDICT: //')
        sv=$(grep -A1 'singular values' "$f" | tail -1 | sed 's/^ *\[ *//')
        # Keep the EXPONENT.  Truncating these to a fixed width silently turns 2.38e-12 into
        # "2.3839893669", which reads as 2.4 and is the difference between a numerical zero and a
        # real singular value.
        r=$(echo "$sv" | cut -d, -f1 | tr -d ' ' | awk '{printf "%.3e", $1}')
        n=$(echo "$sv" | cut -d, -f2 | tr -d ' ' | awk '{printf "%.3e", $1}')
        g=$(grep 'largest ratio jump' "$f" \
            | sed 's/largest ratio jump \([0-9.E+-]*\) after \([0-9]*\).*/\1 (n=\2)/')
    elif grep -qiE 'Runtime error|error,|Segmentation|Killed' "$f" 2>/dev/null; then
        state=ERROR; v=$(grep -iE 'Runtime error|error,' "$f" | head -1 | cut -c1-40); r=; n=; g=
    else
        state=running; v=; r=; n=; g=
        [ ! -s "$f" ] && state='running?'   # 0 bytes: buffered, or killed -- check the PID
    fi
    printf '%-22s %-9s %-18s %-11s %-11s %s\n' "$b" "$state" "$v" "$r" "$n" "$g"
done
echo
echo "magma stage-3 jobs still alive: $(ps -eo command | grep -c '[m]agma.exe -b DISC')"
