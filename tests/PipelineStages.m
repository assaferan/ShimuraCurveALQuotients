// The stage order in workingcode.m (FILTER_STAGES) and run_pipeline.sh: the new stages sit where
// they were placed, and the existing stage names (which are data file names) are unchanged.
// Text-level check: FILTER_STAGES is package-local, so it is read from the source.

function StageOrder(text, names)
    pos := [Position(text, "\"" cat n cat "\"") : n in names];
    return pos, &and[p gt 0 : p in pos] and &and[pos[i] lt pos[i+1] : i in [1..#pos-1]];
end function;

wc := [
    "FilterByTraceStar", "HHProposition1", "FilterByTwistedTraceStar",
    "SpecialFiberIsomorphismStar", "FilterByWeilPolynomialStar",
    "FilterByTwistedWeilPolynomialStar", "FilterStarCurvesByFpAutomorphisms",
    "FilterByNonALInvolutionsStar", "UpdateByGenus", "UpdateCurves1",
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
    "FilterByTraceStar", "HHProposition1", "FilterByTwistedTraceStar",
    "SpecialFiberIsomorphismStar", "FilterByWeilPolynomialStar", "FilterByTwistedWeilPolynomialStar", "FilterStarCurvesByFpAutomorphisms",
    "FilterByNonALInvolutionsStar", "UpdateCurves5", "FilterByAutomorphismGroup",
    "UpdateCurvesAfterAutomorphismGroup", "FilterByTrace", "UpdateCurves6", "FilterByTwistedTrace",
    "UpdateCurvesAfterTwistedTrace", "FilterByWeilPolynomial", "UpdateCurves7",
    "FilterByTwistedWeilPolynomial", "UpdateCurvesAfterTwistedWeilPolynomial",
    "FilterByNonALInvolutions", "UpdateCurves8"];
pos, ok := StageOrder(body, sh);
printf "  run_pipeline.sh positions %o\n", pos;
assert ok;

// The whole star phase, stage by stage, must be the SAME list in both: attribution is the first
// stage to decide a curve, so a star stage that runs at a different point in the two orders labels
// curves differently, and a sequential recompute_data would then fail its equality check against
// data made by run_pipeline.sh.  (The star twisted Weil stage decides 144, 152, 160 and 312, which
// FilterByWeilPolynomialStar decides first; that stage used to be missing from FILTER_STAGES.)
// workingcode.m: the uncommented <"Name", ...> entries of FILTER_STAGES before "UpdateByGenus".
// run_pipeline.sh: the run_seq/run_par stages before the expansion, less FindPairs, which
// compute_data runs outside the stage list.
function Uncommented(line)
    c := Position(line, "//");
    return c eq 0 select line else line[1..c-1];
end function;
wct := Read("workingcode.m");
wct := wct[Position(wct, "FILTER_STAGES := [*")..#wct];
wstar := [];
for line in Split(wct, "\n") do
    ok, _, sub := Regexp("^ *<\"([A-Za-z0-9_]+)\",", Uncommented(line));
    if not ok then continue; end if;
    if sub[1] eq "UpdateByGenus" then break; end if;
    Append(~wstar, sub[1]);
end for;
sstar := [];
for line in Split(body, "\n") do
    ok, _, sub := Regexp("^run_(seq|par) +\"([A-Za-z0-9_]+)\"", line);
    if not ok then continue; end if;
    if sub[2] eq "GetQuotientsAndGenera_UpdateByGenus" then break; end if;
    if sub[2] ne "FindPairs" then Append(~sstar, sub[2]); end if;
end for;
printf "  star phase, workingcode.m:    %o\n", wstar;
printf "  star phase, run_pipeline.sh: %o\n", sstar;
assert #wstar eq 10;
assert wstar eq sstar;

// reconstruct_attribution.m credits each curve to the FIRST stage that decided it, so its
// star_stages cat full_stages must be exactly the uncommented FILTER_STAGES less UpdateGenera and
// the check-only stages, which decide nothing.  A missing stage credits its curves to the next
// UpdateCurves snapshot.  A check-only stage must still RUN (be in FILTER_STAGES and the star
// phase above), and must run its Check intrinsic, which leaves the curves unchanged.
check_only := ["HHProposition1"];
wall := [];
for line in Split(wct[1..Position(wct, "*];")], "\n") do
    ok, _, sub := Regexp("^ *<\"([A-Za-z0-9_]+)\", *([A-Za-z0-9_]+)", Uncommented(line));
    if not ok or sub[1] eq "UpdateGenera" then continue; end if;
    if sub[1] in check_only then
        assert sub[2] eq "Check" cat sub[1];
    else
        Append(~wall, sub[1]);
    end if;
end for;
assert &and[s in wstar : s in check_only];
ra := Read("reconstruct_attribution.m");
ra := ra[Position(ra, "star_stages := [")..Position(ra, "final := Load")];
rall := [];
for line in Split(ra, "\n") do
    ok, _, sub := Regexp("^ *\"([A-Za-z0-9_]+)\",?$", line);
    if ok then Append(~rall, sub[1]); end if;
end for;
printf "  reconstruct_attribution.m: %o stages, FILTER_STAGES: %o\n", #rall, #wall;
assert rall eq wall;

// analysis_stages.m (the per-stage counts table) lists every snapshot that can change a count:
// FILTER_STAGES less the check-only stages, plus FindPairs.
as := Read("analysis_stages.m");
as := as[Position(as, "stages := [")..Position(as, "];")];
asl := [];
for line in Split(as, "\n") do
    ok, _, sub := Regexp("^ *\"([A-Za-z0-9_]+)\",?$", line);
    if ok then Append(~asl, sub[1]); end if;
end for;
assert asl eq ["FindPairs", "UpdateGenera"] cat wall;
// ... and the paper tables give a check-only stage no row.
assert &and[Position(Read("make_latex_tables.m"), "curves_after_" cat s cat ".dat") eq 0 : s in check_only];

// The check-only stage on the committed FilterByTraceStar snapshot: [HH] Table 2 and Proposition 1
// are verified (on the input, and on HHProposition1 run on a copy), and the curves are unchanged.
// ShimuraQuot has reference semantics, so a shallow copy inside CheckHHProposition1 would write
// HH's four verdicts into the input; its own final assert and the comparison below catch that.
star := eval Read("data/curves_after_FilterByTraceStar.dat");
before := Sprint(star, "Magma");
CheckHHProposition1(~star);
assert Sprint(star, "Magma") eq before;
copy := eval before;
HHProposition1(~copy);
hh := [i : i in [1..#copy] | assigned copy[i]`IsSubhyp and not assigned star[i]`IsSubhyp];
assert [<copy[i]`D, copy[i]`N> : i in hh] eq [<1, 194>, <1, 546>, <205, 3>, <1995, 2>];
