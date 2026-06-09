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
NUM_WORKERS="${2:-8}"
DATA_DIR="${3:-data}"
NUM_CHUNKS="${4:-${NUM_WORKERS}}"

MAGMA_CMD="${MAGMA_CMD:-magma}"

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${REPO_DIR}"

CHUNKS_DIR="${DATA_DIR}/parallel_chunks_${STAGE}"

# Map each stage to the .dat file produced by the preceding pipeline step
case "${STAGE}" in
    FilterByTraceStar)
        INPUT_DAT="${DATA_DIR}/curves_after_UpdateByGenusStar.dat"
        ;;
    FilterStarCurvesByFpAutomorphisms)
        INPUT_DAT="${DATA_DIR}/curves_after_HHProposition1.dat"
        ;;
    FilterByALFixedPointsOnQuotient)
        INPUT_DAT="${DATA_DIR}/curves_after_UpdateCurves1.dat"
        ;;
    FilterByDegeneracyMorphism)
        INPUT_DAT="${DATA_DIR}/curves_after_UpdateCurves3.dat"
        ;;
    FilterByComplicatedALFixedPointsOnQuotient)
        INPUT_DAT="${DATA_DIR}/curves_after_UpdateCurves4.dat"
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
        echo "Supported: FilterByTraceStar, FilterStarCurvesByFpAutomorphisms," >&2
        echo "           FilterByALFixedPointsOnQuotient, FilterByDegeneracyMorphism," >&2
        echo "           FilterByComplicatedALFixedPointsOnQuotient, FilterByTrace," >&2
        echo "           FilterByWeilPolynomial, FilterByNonALInvolutions" >&2
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

echo "All ${NUM_CHUNKS} chunk files present. Merging..."

python3 parallel_merge.py "${CHUNKS_DIR}" "${NUM_CHUNKS}" "${OUTPUT_DAT}"

T_END=$(date +%s)
ELAPSED=$((T_END - T_START))
echo ""
echo "Done in ${ELAPSED}s. Output: ${OUTPUT_DAT}"
