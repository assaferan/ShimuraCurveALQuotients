#!/usr/bin/env bash
# Per-group timeout runner for the ratpts sweeps (design 2026-06-14, Section 3).
#
# Runs each requested TABLE index as its OWN `magma idx:=i <driver>` process under
# a wall-clock timeout, so a single hang/OOM (e.g. 10,61) cannot wedge the whole
# sweep. There is no coreutils `timeout`/`gtimeout` on this box, so we use a
# portable perl-`alarm` wrapper that puts magma in its own process group and
# SIGKILLs the whole group (magma + any polymake children) when the alarm fires.
#
# Usage:
#   ./run_table.sh <timeout-seconds> <idx> [idx ...]
#   DRIVER=ratpts_table1.m ./run_table.sh 900 1 2 3 4 24 26 27 28 31 35
#
# Verdict rows (group header + the driver's verdict lines, plus TIMEOUT) are
# appended to $RESULTS; full output is tee'd to a fresh per-run log.
set -u

DRIVER="${DRIVER:-ratpts_table1.m}"
RESULTS="${RESULTS:-ratpts_table1_output.txt}"

if [ "$#" -lt 2 ]; then
    echo "usage: $0 <timeout-seconds> <idx> [idx ...]" >&2
    exit 2
fi
TIMEOUT="$1"; shift
INDICES=("$@")

STAMP=$(date +%Y%m%d_%H%M%S)
LOG="table1plus7_plausible_${STAMP}.log"

banner="=== plausible sweep: driver=$DRIVER timeout=${TIMEOUT}s indices=[${INDICES[*]}] start=$(date) ==="
echo "$banner" | tee -a "$LOG"
echo "$banner" >> "$RESULTS"

for i in "${INDICES[@]}"; do
    tmp=$(mktemp)
    header="########## idx=$i  $(date '+%Y-%m-%d %H:%M:%S') ##########"
    echo "" | tee -a "$LOG"
    echo "$header" | tee -a "$LOG"
    start=$(date +%s)

    # perl-alarm wrapper: fork, child becomes its own process-group leader and
    # execs magma; parent SIGKILLs the whole group on alarm and exits 124.
    perl -e '
        my $t = shift;
        my $pid = fork();
        if ($pid == 0) { setpgrp(0,0); exec @ARGV or die "exec failed: $!\n"; }
        $SIG{ALRM} = sub { kill("-KILL", $pid); waitpid($pid, 0); exit 124; };
        alarm $t;
        waitpid($pid, 0);
        exit($? >> 8);
    ' "$TIMEOUT" magma idx:="$i" "$DRIVER" 2>&1 | tee -a "$LOG" "$tmp"
    rc=${PIPESTATUS[0]}
    dur=$(( $(date +%s) - start ))

    {
        echo ""
        echo "$header (dur=${dur}s exit=$rc)"
        grep -E '==== D=|model y\^2|IS elliptic|E \(min model\)|aInvariants|genus 1, IsEllipticCurve found NO|insufficient CM|not an immediate cover|not among computed|WARNING genus|ERROR on|computed .* cover' "$tmp"
        if [ "$rc" -eq 124 ]; then
            echo "  TIMEOUT(>${TIMEOUT}s) -- killed, recorded, not retried"
        fi
    } >> "$RESULTS"

    if [ "$rc" -eq 124 ]; then
        echo "  >>> idx=$i TIMEOUT after ${dur}s (>${TIMEOUT}s) -- killed group, moving on" | tee -a "$LOG"
    else
        echo "  >>> idx=$i exit=$rc dur=${dur}s" | tee -a "$LOG"
    fi
    rm -f "$tmp"
done

echo "" | tee -a "$LOG"
echo "=== sweep done $(date) ===" | tee -a "$LOG"
echo "=== sweep done $(date) ===" >> "$RESULTS"
