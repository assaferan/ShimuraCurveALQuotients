// Recompute per-stage open/ruled/proved counts across the pipeline output files.
// Writes a table to data/stage_counts.txt.
SetQuitOnError(true);
AttachSpec("ShimuraQuotients.spec");
SetVerbose("ShimuraQuotients", 0);

out := "data/stage_counts.txt";
procedure W(s) Write(out, s); end procedure;

function Status(X)
    if not assigned X`IsSubhyp then return "open"; end if;
    if X`IsSubhyp then return "proved"; else return "ruled"; end if;
end function;

stages := [
    "FindPairs",
    "UpdateGenera",
    "UpdateByGenusStar",
    "FilterByTraceStar",
    "FilterByTwistedTraceStar",
    "HHProposition1",
    "SpecialFiberIsomorphismStar",
    "FilterByWeilPolynomialStar",
    "FilterByTwistedWeilPolynomialStar",
    "FilterStarCurvesByFpAutomorphisms",
    "FilterByNonALInvolutionsStar",
    "UpdateByGenus",
    "UpdateCurves1",
    "FilterByALFixedPointsOnQuotient",
    "UpdateCurves2",
    "Genus3CoversGenus2",
    "UpdateCurves3",
    "FilterByDegeneracyMorphism",
    "UpdateCurves4",
    "FilterByComplicatedALFixedPointsOnQuotient",
    "FilterByGeneralizedComplicatedFixedPoints",
    "UpdateCurves5",
    "FilterByAutomorphismGroup",
    "UpdateCurvesAfterAutomorphismGroup",
    "FilterByTrace",
    "UpdateCurves6",
    "FilterByTwistedTrace",
    "UpdateCurvesAfterTwistedTrace",
    "FilterByWeilPolynomial",
    "UpdateCurves7",
    "FilterByTwistedWeilPolynomial",
    "UpdateCurvesAfterTwistedWeilPolynomial",
    "FilterByNonALInvolutions",
    "UpdateCurves8"
];

W(Sprintf("%-3o %-44o %-6o %-6o %-7o %-6o %-6o", "#", "Step", "total", "open", "ruled", "proved", "+ruled"));
// A stage whose file does not exist yet (the stages added after the last full run, until the
// rerun) is left out of the table, and the remaining rows are numbered consecutively; on data
// from before those stages the table is therefore unchanged.
function Exists(f)
    try
        _ := Open(f, "r");
        return true;
    catch e
        return false;
    end try;
end function;

prevR := 0;
i := 0;
for s in stages do
    f := "data/curves_after_" cat s cat ".dat";
    if not Exists(f) then
        printf "skipping %o: %o does not exist (stage not run yet)\n", s, f;
        continue;
    end if;
    i +:= 1;
    curves := eval Read(f);
    n := #curves;
    o := #[X : X in curves | Status(X) eq "open"];
    r := #[X : X in curves | Status(X) eq "ruled"];
    p := #[X : X in curves | Status(X) eq "proved"];
    W(Sprintf("%-3o %-44o %-6o %-6o %-7o %-6o %-6o", i, s, n, o, r, p, r-prevR));
    prevR := r;
    delete curves;
end for;

printf "wrote %o\n", out;
quit;
