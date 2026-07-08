#!/usr/bin/env bash
# Parallelises a single FilterBy* stage across multiple Magma processes
# using GNU parallel.
#
# Usage:
#   ./run_parallel_filter.sh <stage> [num_workers] [data_dir] [num_chunks]
#
# num_workers  - number of Magma processes to run simultaneously (default: 8)
# num_chunks   - number of pieces to split the curve list into (default: num_workers)
#                Set num_chunks > num_workers to allow finer-grained checkpointing.
# data_dir     - defaults to "data".  The output is written to
#                <data_dir>/curves_after_<stage>.dat

set -euo pipefail

STAGE="${1:-FilterByTrace}"
NUM_WORKERS="${2:-128}"
DATA_DIR="${3:-data}"
# Many more chunks than workers: with the cost-aware worker, fine chunks isolate the heavy
# curves so the makespan approaches the slowest single curve rather than the slowest cluster.
NUM_CHUNKS="${4:-1024}"

MAGMA_CMD="${MAGMA_CMD:-magma}"

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${REPO_DIR}"

CHUNKS_DIR="${DATA_DIR}/parallel_chunks_${STAGE}"

# Map each stage to the .dat file produced by the preceding pipeline step
case "${STAGE}" in
    FilterByTraceStar)
        INPUT_DAT="${DATA_DIR}/curves_after_UpdateByGenusStar.dat"
        ;;
    FilterByWeilPolynomialStar)
        # Weil-polynomial filter run on the star curves, before FilterStarCurvesByFpAutomorphisms.
        # The Fp-automorphism ramification count and the Weil polynomial both derive from the same
        # F_{p^d} point counts, and the Weil filter is at least as strong (verified: it rules out
        # every star curve FpAut does), so this stage subsumes the FpAut stage that follows.
        INPUT_DAT="${DATA_DIR}/curves_after_SpecialFiberIsomorphismStar.dat"
        ;;
    FilterStarCurvesByFpAutomorphisms)
        INPUT_DAT="${DATA_DIR}/curves_after_FilterByWeilPolynomialStar.dat"
        ;;
    FilterByNonALInvolutionsStar)
        # Non-AL involution filter run on the star curves (full AL group), before they
        # are expanded into all quotients.  A star curve proven non-subhyperelliptic
        # propagates to its whole cover-tree via UpwardClosure after expansion.
        INPUT_DAT="${DATA_DIR}/curves_after_FilterStarCurvesByFpAutomorphisms.dat"
        ;;
    FilterByALFixedPointsOnQuotient)
        # Preceded by FilterBySpecialFiber (moved to the front of the all-quotients block so its
        # cheap determinations prune every downstream filter), hence reads its output.
        INPUT_DAT="${DATA_DIR}/curves_after_FilterBySpecialFiber.dat"
        ;;
    FilterByDegeneracyMorphism)
        INPUT_DAT="${DATA_DIR}/curves_after_UpdateCurves3.dat"
        ;;
    FilterByComplicatedALFixedPointsOnQuotient)
        INPUT_DAT="${DATA_DIR}/curves_after_UpdateCurves4.dat"
        ;;
    FilterByGeneralizedComplicatedFixedPoints)
        # Generalized [FH] Prop 6 via mixed groups <W_odd, V_p>; additive to the AL version above.
        INPUT_DAT="${DATA_DIR}/curves_after_FilterByComplicatedALFixedPointsOnQuotient.dat"
        ;;
    FilterBySpecialFiber)
        # Special-fiber (reduction mod p) non-hyperellipticity test, [FH] Section 5 generalized.
        # D=1, D=6, D=10 and D=22 curves with a genus-0 special-fiber component.  Runs first in the
        # all-quotients block (right after UpdateCurves1): it is by far the cheapest filter, so
        # its determinations prune all the heavier filters that follow.
        INPUT_DAT="${DATA_DIR}/curves_after_UpdateCurves1.dat"
        ;;
    FilterByTrace)
        INPUT_DAT="${DATA_DIR}/curves_after_UpdateCurves5.dat"
        ;;
    FilterByWeilPolynomial)
        INPUT_DAT="${DATA_DIR}/curves_after_UpdateCurves6.dat"
        ;;
    FilterByNonALInvolutions)
        INPUT_DAT="${DATA_DIR}/curves_after_UpdateCurves7.dat"
        ;;
    *)
        echo "ERROR: unknown stage '${STAGE}'" >&2
        echo "Supported: FilterByTraceStar, FilterByWeilPolynomialStar," >&2
        echo "           FilterStarCurvesByFpAutomorphisms," >&2
        echo "           FilterByALFixedPointsOnQuotient, FilterByDegeneracyMorphism," >&2
        echo "           FilterByComplicatedALFixedPointsOnQuotient," >&2
        echo "           FilterByGeneralizedComplicatedFixedPoints, FilterBySpecialFiber," >&2
        echo "           FilterByTrace," >&2
        echo "           FilterByWeilPolynomial, FilterByNonALInvolutions," >&2
        echo "           FilterByNonALInvolutionsStar" >&2
        exit 1
        ;;
esac

OUTPUT_DAT="${DATA_DIR}/curves_after_${STAGE}.dat"

if [ ! -f "${INPUT_DAT}" ]; then
    echo "ERROR: input file not found: ${INPUT_DAT}" >&2
    exit 1
fi

JOBLOG="${CHUNKS_DIR}/joblog.txt"
RESULTS_DIR="${CHUNKS_DIR}/results"

mkdir -p "${CHUNKS_DIR}" "${RESULTS_DIR}"

echo "========================================"
echo "Stage:    ${STAGE}"
echo "Workers:  ${NUM_WORKERS} (parallel jobs)"
echo "Chunks:   ${NUM_CHUNKS} (data splits)"
echo "Input:    ${INPUT_DAT}"
echo "Output:   ${OUTPUT_DAT}"
echo "ChunkDir: ${CHUNKS_DIR}"
echo "JobLog:   ${JOBLOG}"
echo "Results:  ${RESULTS_DIR}/ (per-chunk stdout/stderr)"
echo "========================================"

T_START=$(date +%s)

parallel \
    --jobs "${NUM_WORKERS}" \
    --joblog "${JOBLOG}" \
    --results "${RESULTS_DIR}/chunk_{}" \
    --halt now,fail=1 \
    --eta \
    ${MAGMA_CMD} \
        "input_dat:=${INPUT_DAT}" \
        "chunk:={}" \
        "total_chunks:=${NUM_CHUNKS}" \
        "stage:=${STAGE}" \
        "output_dat:=${CHUNKS_DIR}/chunk_{}_of_${NUM_CHUNKS}.dat" \
        parallel_filter_worker.m \
    ::: $(seq 1 "${NUM_CHUNKS}")

echo ""
echo "All workers finished. Verifying chunk files..."

MISSING=0
for c in $(seq 1 "${NUM_CHUNKS}"); do
    CHUNK_FILE="${CHUNKS_DIR}/chunk_${c}_of_${NUM_CHUNKS}.dat"
    if [ ! -f "${CHUNK_FILE}" ]; then
        echo "  ERROR: missing chunk file: ${CHUNK_FILE}" >&2
        MISSING=$((MISSING + 1))
    fi
done
if [ "${MISSING}" -gt 0 ]; then
    echo "ERROR: ${MISSING} chunk file(s) missing; aborting merge." >&2
    exit 1
fi

echo "All ${NUM_CHUNKS} chunk files present. Merging (by original index)..."

# The chunk files are index-tagged, so the merge reorders by original index (done in Magma).
${MAGMA_CMD} \
    "chunks_dir:=${CHUNKS_DIR}" \
    "total_chunks:=${NUM_CHUNKS}" \
    "output_dat:=${OUTPUT_DAT}" \
    parallel_merge.m

T_END=$(date +%s)
ELAPSED=$((T_END - T_START))
echo ""
echo "Done in ${ELAPSED}s. Output: ${OUTPUT_DAT}"
