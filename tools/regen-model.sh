#!/bin/bash
# tools/regen-model.sh -- regenerate a committed model with the flags it actually needs.
#
# WHY THIS EXISTS. Five committed models do not regenerate under default settings, for two
# different reasons, and before 2026-09-05 nothing recorded what they DID need. Running
# genmodels.m by hand and getting a smaller model back looks like a regression and is not one.
# data/models/PROVENANCE.md is the source of truth; this script is that table made executable.
#
#   ./tools/regen-model.sh 39 2 /tmp/out      # picks up CMNONCOPRIME=1 automatically
#   ./tools/regen-model.sh 51 1 /tmp/out      # no flags needed
#
# Then diff against the committed file (ignoring the header comment):
#   diff <(tail -n +2 data/models/models_39_2.m) <(tail -n +2 /tmp/out/models_39_2.m)
set -u

D="${1:?usage: regen-model.sh D N OUTDIR}"
N="${2:?usage: regen-model.sh D N OUTDIR}"
OUT="${3:?usage: regen-model.sh D N OUTDIR}"

: "${NORMALIZ_BIN:=$HOME/Documents/GitHub/normaliz-3.11.1/normaliz}"
export NORMALIZ_BIN
if [ ! -x "$NORMALIZ_BIN" ]; then
    echo "ERROR: NORMALIZ_BIN=$NORMALIZ_BIN is not executable." >&2
    echo "  Without it a polytope solve fails SILENTLY -- you get 'no solutions', not an error." >&2
    exit 1
fi

# genmodels.m is triage tooling and lives on the campaign branch, not on main (see CLAUDE.md).
GM=""
for cand in "worktrees/campaign/vvdata/weyl-campaign/genmodels.m" \
            "vvdata/weyl-campaign/genmodels.m" \
            "genmodels.m"; do
    [ -f "$cand" ] && { GM="$cand"; break; }
done
if [ -z "$GM" ]; then
    echo "ERROR: genmodels.m not found. It lives on the m0-theta-campaign branch at" >&2
    echo "  vvdata/weyl-campaign/genmodels.m -- add a worktree, e.g." >&2
    echo "  git worktree add worktrees/campaign m0-theta-campaign" >&2
    exit 1
fi

# --- the flags each base needs, mirroring data/models/PROVENANCE.md -------------------------
# Keep these two in sync. A base missing from this table regenerates with no flags.
FLAGS=""
case "${D}_${N}" in
    39_2|14_3)      FLAGS="CMNONCOPRIME=1" ;;   # deliberate: no theoretical guarantee, oracle-validated
    22_5|15_2|22_3) FLAGS="Y2TWIST=1" ;;        # accidental: the y2-scale guard postdates these files
esac

mkdir -p "$OUT"
echo "regenerating X_0^${D}(${N})${FLAGS:+ with $FLAGS}"
[ -n "$FLAGS" ] && echo "  (see data/models/PROVENANCE.md for why)"

# BFPROGRESS/COVPROGRESS are unbuffered (WriteStderr), so a killed run keeps its progress.
# Harmless when the base is fast; essential on the ones that take hours.
env $FLAGS BFPROGRESS=1 COVPROGRESS=1 \
    magma -b "D_s:=${D}" "N_s:=${N}" "OUTDIR:=${OUT}" "$GM" < /dev/null \
    > "${OUT}/${D}_${N}.log" 2>&1
rc=$?

F="${OUT}/models_${D}_${N}.m"
if [ -f "$F" ]; then
    echo "  wrote $F"
    C="data/models/models_${D}_${N}.m"
    if [ -f "$C" ]; then
        # Compare the models[...] lines only. Committed files carry a provenance header of
        # varying length, so "tail -n +2" does NOT align them -- that mis-comparison reported a
        # false DIFFERS the first time this script ran.
        if diff -q <(grep '^models\[' "$C" | sort) <(grep '^models\[' "$F" | sort) >/dev/null 2>&1; then
            echo "  IDENTICAL to the committed model (cover lines match)"
        else
            echo "  DIFFERS from the committed model:"
            diff <(grep '^models\[' "$C" | sort) <(grep '^models\[' "$F" | sort) | head -20 | sed 's/^/    /'
            echo "  (a re-presented but isomorphic curve is not a regression; a MISSING cover is)"
        fi
    fi
else
    echo "  NO MODEL WRITTEN (exit $rc). Last lines:"
    tail -5 "${OUT}/${D}_${N}.log" | sed 's/^/    /'
fi
exit $rc
