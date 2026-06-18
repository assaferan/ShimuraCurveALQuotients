#!/usr/bin/env bash
# Build data/pipeline_summary.txt and its census/open analyses from the
# completed-run data files (data/curves_after_*.dat).
#
# pipeline_summary.txt has no single generator: it is the hand-maintained
# front-matter (data/pipeline_summary_frontmatter.txt -- algorithm description,
# Phase A/B stage tables, expansion breakdown, END OF PIPELINE totals, attribution
# notes) concatenated with two machine-generated census reports.
#
# The stage tables inside the front-matter are pasted by hand from
# data/stage_counts.txt (regenerated below); update them there if the data changes.
#
# Usage:  ./make_pipeline_summary.sh
set -euo pipefail
cd "$(dirname "$0")"

FRONT=data/pipeline_summary_frontmatter.txt
CENSUS=data/pipeline_analysis_census.txt
OPEN=data/open_cases_analysis.txt
STAGES=data/stage_counts.txt
OUT=data/pipeline_summary.txt

# 1. per-stage open/ruled/proved counts (paste these into the front-matter tables)
rm -f "$STAGES"
magma analysis_stages.m

# 2. machine-generated census + open-case analyses (both Write-append, so wipe first)
rm -f "$CENSUS"
magma analysis_census.m
rm -f "$OPEN"
magma analysis_open.m

# 3. assemble: hand front-matter (ends at the "DETAILED CENSUS" header) + census + open
cat "$FRONT" "$CENSUS" "$OPEN" > "$OUT"
echo "wrote $OUT  ($(wc -l < "$OUT") lines)"
