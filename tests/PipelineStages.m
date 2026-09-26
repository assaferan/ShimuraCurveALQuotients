// The stage order in workingcode.m (FILTER_STAGES) and run_pipeline.sh: the new stages sit where
// they were placed, and the existing stage names (which are data file names) are unchanged.
// Text-level check: FILTER_STAGES is package-local, so it is read from the source.

function StageOrder(text, names)
    pos := [Position(text, "\"" cat n cat "\"") : n in names];
    return pos, &and[p gt 0 : p in pos] and &and[pos[i] lt pos[i+1] : i in [1..#pos-1]];
end function;

wc := [
    "FilterByTraceStar", "FilterByTwistedTraceStar", "FilterByTwistedWeilPolynomialStar",
    "HHProposition1", "FilterByNonALInvolutionsStar", "UpdateByGenus", "UpdateCurves1",
    "FilterByComplicatedALFixedPointsOnQuotient", "FilterByGeneralizedComplicatedFixedPoints",
    "UpdateCurves5", "FilterByAutomorphismGroup", "UpdateCurvesAfterAutomorphismGroup",
    "FilterByTrace", "UpdateCurves6", "FilterByTwistedTrace", "UpdateCurvesAfterTwistedTrace",
    "FilterByWeilPolynomial", "UpdateCurves7", "FilterByTwistedWeilPolynomial",
    "UpdateCurvesAfterTwistedWeilPolynomial", "FilterByNonALInvolutions", "UpdateCurves8"];
pos, ok := StageOrder(Read("workingcode.m"), wc);
printf "  workingcode.m positions %o\n", pos;
assert ok;
// the final stage is still UpdateCurves8, so GetHyperellipticCandidates reads the committed data
assert Position(Read("workingcode.m"), "<\"UpdateCurves8\", UpdateCurves>\n*];") gt 0;

rp := Read("run_pipeline.sh");
body := rp[Position(rp, "set -euo pipefail")..#rp];   // skip the header comment
sh := [
    "FilterByTraceStar", "FilterByTwistedTraceStar", "HHProposition1", "FilterByWeilPolynomialStar",
    "FilterByTwistedWeilPolynomialStar", "FilterStarCurvesByFpAutomorphisms",
    "FilterByNonALInvolutionsStar", "UpdateCurves5", "FilterByAutomorphismGroup",
    "UpdateCurvesAfterAutomorphismGroup", "FilterByTrace", "UpdateCurves6", "FilterByTwistedTrace",
    "UpdateCurvesAfterTwistedTrace", "FilterByWeilPolynomial", "UpdateCurves7",
    "FilterByTwistedWeilPolynomial", "UpdateCurvesAfterTwistedWeilPolynomial",
    "FilterByNonALInvolutions", "UpdateCurves8"];
pos, ok := StageOrder(body, sh);
printf "  run_pipeline.sh positions %o\n", pos;
assert ok;
