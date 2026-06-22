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
    "HHProposition1",
    "SpecialFiberIsomorphismStar",
    "FilterByWeilPolynomialStar",
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
    "FilterByTrace",
    "UpdateCurves6",
    "FilterByWeilPolynomial",
    "UpdateCurves7",
    "FilterByNonALInvolutions",
    "UpdateCurves8"
];

W(Sprintf("%-3o %-44o %-6o %-6o %-7o %-6o %-6o", "#", "Step", "total", "open", "ruled", "proved", "+ruled"));
prevR := 0;
for i in [1..#stages] do
    s := stages[i];
    f := "data/curves_after_" cat s cat ".dat";
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
