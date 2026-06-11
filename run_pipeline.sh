#!/usr/bin/env bash
# Full parallelized pipeline — equivalent to workingcode.m's GetHyperellipticCandidates
# but with all FilterBy* stages run in parallel across NUM_WORKERS Magma processes.
#
# Usage:
#   ./run_pipeline.sh [num_workers] [data_dir] [num_chunks]
#
# num_workers - parallel Magma jobs at once (default: 8)
# num_chunks  - pieces to split the curve list into (default: num_workers)
# data_dir    - defaults to "data/par".  All intermediate and final .dat files
#               are written there, leaving the sequential pipeline's "data/"
#               directory untouched.
#
# Stages and their parallelisation:
#
#   FindPairs                              [sequential]
#   UpdateGenera                           [sequential]  + VerifyHHTable1
#   UpdateByGenusStar                      [sequential]
#   FilterByTraceStar                      [PARALLEL]
#   HHProposition1                         [sequential]  + VerifyHHTable2 (pre), VerifyHHProposition1 (post)
#   FilterStarCurvesByFpAutomorphisms      [PARALLEL]
#   FilterByNonALInvolutionsStar           [PARALLEL]  (star curves, pre-expansion)
#   GetQuotientsAndGenera + UpdateByGenus  [sequential]  (carries star determinations
#                                                         onto the full-W entries)
#                                                         + VerifyFHTheorem3 (post)
#   UpdateCurves1                          [sequential]
#   FilterByALFixedPointsOnQuotient        [PARALLEL]
#   UpdateCurves2                          [sequential]
#   Genus3CoversGenus2                     [sequential]
#   UpdateCurves3                          [sequential]
#   FilterByDegeneracyMorphism             [PARALLEL]
#   UpdateCurves4                          [sequential]
#   FilterByComplicatedALFixedPointsOnQuotient [PARALLEL]
#   UpdateCurves5                          [sequential]  + VerifyFHTable3 (post)
#   FilterByTrace                          [PARALLEL]
#   UpdateCurves6                          [sequential]
#   FilterByWeilPolynomial                 [PARALLEL]
#   UpdateCurves7                          [sequential]
#   FilterByNonALInvolutions               [PARALLEL]
#   UpdateCurves8                          [sequential]

set -euo pipefail

# 128 workers per policy on this shared machine; many fine chunks so the cost-aware worker
# can isolate the heavy tail curves (makespan -> slowest single curve, not slowest cluster).
NUM_WORKERS="${1:-128}"
DATA_DIR="${2:-data/par}"
NUM_CHUNKS="${3:-1024}"
MAGMA_CMD="${MAGMA_CMD:-magma -b}"
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${REPO_DIR}"

mkdir -p "${DATA_DIR}"

# ── helpers ──────────────────────────────────────────────────────────────────

run_seq() {
    local stage="$1"
    local input_dat="${2:-}"
    local output_dat="$3"

    if [ -f "${output_dat}" ]; then
        echo "[skip] ${stage} — ${output_dat} already exists"
        return
    fi

    echo ""
    echo "=== [sequential] ${stage} ==="
    if [ -n "${input_dat}" ]; then
        ${MAGMA_CMD} \
            "stage:=${stage}" \
            "input_dat:=${input_dat}" \
            "output_dat:=${output_dat}" \
            run_sequential_stage.m
    else
        ${MAGMA_CMD} \
            "stage:=${stage}" \
            "output_dat:=${output_dat}" \
            run_sequential_stage.m
    fi
}

run_par() {
    local stage="$1"
    local output_dat="${DATA_DIR}/curves_after_${stage}.dat"

    if [ -f "${output_dat}" ]; then
        echo "[skip] ${stage} — ${output_dat} already exists"
        return
    fi

    echo ""
    echo "=== [parallel jobs=${NUM_WORKERS} chunks=${NUM_CHUNKS}] ${stage} ==="
    MAGMA_CMD="${MAGMA_CMD}" ./run_parallel_filter.sh "${stage}" "${NUM_WORKERS}" "${DATA_DIR}" "${NUM_CHUNKS}"
}

# ── star-curve stages (2342 curves) ──────────────────────────────────────────

D="${DATA_DIR}"

run_seq "FindPairs"           ""                                                   "${D}/curves_after_FindPairs.dat"
run_seq "UpdateGenera"        "${D}/curves_after_FindPairs.dat"                    "${D}/curves_after_UpdateGenera.dat"
run_seq "UpdateByGenusStar"   "${D}/curves_after_UpdateGenera.dat"                 "${D}/curves_after_UpdateByGenusStar.dat"
run_par "FilterByTraceStar"
run_seq "HHProposition1"      "${D}/curves_after_FilterByTraceStar.dat"            "${D}/curves_after_HHProposition1.dat"
run_par "FilterStarCurvesByFpAutomorphisms"
# Non-AL involution filter on the star curves themselves: a star curve proven
# non-subhyperelliptic here prunes its entire cover-tree (via UpwardClosure) right
# after expansion, sparing all its quotients the heavy per-curve stages below.
run_par "FilterByNonALInvolutionsStar"

# ── expand to all AL quotients ────────────────────────────────────────────────

run_seq "GetQuotientsAndGenera_UpdateByGenus" \
                              "${D}/curves_after_FilterByNonALInvolutionsStar.dat" \
                              "${D}/curves_after_UpdateByGenus.dat"

# ── all-quotients stages (~18 000+ curves) ────────────────────────────────────

run_seq "UpdateCurves1"       "${D}/curves_after_UpdateByGenus.dat"                "${D}/curves_after_UpdateCurves1.dat"
run_par "FilterByALFixedPointsOnQuotient"
run_seq "UpdateCurves2"       "${D}/curves_after_FilterByALFixedPointsOnQuotient.dat" \
                                                                                   "${D}/curves_after_UpdateCurves2.dat"
run_seq "Genus3CoversGenus2"  "${D}/curves_after_UpdateCurves2.dat"                "${D}/curves_after_Genus3CoversGenus2.dat"
run_seq "UpdateCurves3"       "${D}/curves_after_Genus3CoversGenus2.dat"           "${D}/curves_after_UpdateCurves3.dat"
run_par "FilterByDegeneracyMorphism"
run_seq "UpdateCurves4"       "${D}/curves_after_FilterByDegeneracyMorphism.dat"   "${D}/curves_after_UpdateCurves4.dat"
run_par "FilterByComplicatedALFixedPointsOnQuotient"
run_seq "UpdateCurves5"       "${D}/curves_after_FilterByComplicatedALFixedPointsOnQuotient.dat" \
                                                                                   "${D}/curves_after_UpdateCurves5.dat"
run_par "FilterByTrace"
run_seq "UpdateCurves6"       "${D}/curves_after_FilterByTrace.dat"                "${D}/curves_after_UpdateCurves6.dat"
run_par "FilterByWeilPolynomial"
run_seq "UpdateCurves7"       "${D}/curves_after_FilterByWeilPolynomial.dat"       "${D}/curves_after_UpdateCurves7.dat"
run_par "FilterByNonALInvolutions"
run_seq "UpdateCurves8"       "${D}/curves_after_FilterByNonALInvolutions.dat"     "${D}/curves_after_UpdateCurves8.dat"

echo ""
echo "Pipeline complete. Final output: ${D}/curves_after_UpdateCurves8.dat"
