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
// compute_data has no VerifyHHTable2 of its own after FilterByTraceStar: the next stage,
// CheckHHProposition1, runs it on the same curves.
assert wstar[Position(wstar, "FilterByTraceStar") + 1] eq "HHProposition1";

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
// ... and run_pipeline.sh runs it through run_sequential_stage.m, whose case must also call the
// Check intrinsic, not the stage's own (which would write its verdicts into the curves).
rs := Read("run_sequential_stage.m");
// Only code counts: /* */ and // comments are stripped first, and only what is left is tested.
function StripBlockComments(t)
    while Position(t, "/*") gt 0 do
        a := Position(t, "/*");
        b := Position(t[a+2..#t], "*/");
        t := t[1..a-1] cat (b eq 0 select "" else t[a+b+3..#t]);
    end while;
    return t;
end function;
// The whole file is stripped before the case is located, so a comment cannot move its boundaries.
rscode := [Uncommented(l) : l in Split(StripBlockComments(rs), "\n")];
for s in check_only do
    label := "^ *when \"" cat s cat "\": *";
    at := [j : j in [1..#rscode] | Regexp(label, rscode[j])];
    assert #at eq 1;
    // the case body: the rest of the label line, up to the next when / else / end case of this
    // case (indented no deeper than the label, so an if-else inside the body does not end it)
    _, lbl := Regexp(label, rscode[at[1]]);
    code := [rscode[at[1]][#lbl+1..#rscode[at[1]]]];
    ind := Position(lbl, "when") - 1;
    for j in [at[1]+1..#rscode] do
        _, lead := Regexp("^ *", rscode[j]);
        if #lead le ind and Regexp("^ *(when |else|end case)", rscode[j]) then break; end if;
        Append(~code, rscode[j]);
    end for;
    calls := [l : l in code | Regexp("[A-Za-z0-9_]+ *\\( *~ *curves *\\)", l)];
    printf "  run_sequential_stage.m %o: %o\n", s, calls;
    assert #calls eq 1 and Regexp("^ *Check" cat s cat " *\\( *~ *curves *\\) *; *$", calls[1]);
    // ... and no bare call of the stage itself anywhere in the body
    assert not exists{l : l in code | Regexp("(^|[^A-Za-z0-9_])" cat s cat "([^A-Za-z0-9_]|$)", l)};
end for;
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

// Restricted to genus <= 12 (no undecided D = 1 curves above), GetModularByGenus is shorter than
// Table 2 and both checks still pass.  Without HH's marks, 194 and 546 are reported in genus 3, 4.
// 12: the largest genus in GetHHTable2 (N = 2310), which omits [HH] Table 2's genus-19 N = 1680.
// 194 has genus 3 and 546 genus 4 ([HH] Table 2, p. 183); both are removed by Prop. 1, p. 184.
low := [X : X in copy | X`g le 12];
VerifyHHTable2([X : X in star | X`g le 12]);
VerifyHHProposition1(low);
msg := "";
try
    VerifyHHProposition1([X : X in star | X`g le 12]);
catch e
    msg := e`Object;
end try;
assert Position(msg, "in genus [ 3, 4 ]") gt 0;

// An uncertified source is reported only if its target stays unmarked.  Genuine star curves, but the
// targets' genera are hand-set to the sources' so that HHProposition1 applies:
//  * X_0^*(1290): X_0^*(258) at p = 5, not certified, is recorded first; X_0^*(430) at p = 3 then marks it.
//  * X_0^*(1120): only X_0^*(160) at p = 7, not certified, so it stays unmarked and is reported.
// Source genera 3, 3, 4: dim of the all-signs-+1 AL subspace of S_2(Gamma_0(N)), Magma modular symbols.
// Target genera hand-set to 3, 4; the true genera of X_0^*(1290), X_0^*(1120) are 10, 19.
mkstar := function(N, g, id)
    X := CreateShimuraQuot(1, N, {d : d in Divisors(N) | GCD(d, N div d) eq 1});
    X`g := g;
    X`CurveID := id;
    return X;
end function;
two := [mkstar(258, 3, 1), mkstar(430, 3, 2), mkstar(160, 4, 3), mkstar(1290, 3, 4), mkstar(1120, 4, 5)];
for i in [1..3] do
    assert GenusShimuraCurveQuotient(1, two[i]`N, two[i]`W) eq two[i]`g;
    two[i]`IsSubhyp := false; two[i]`IsHyp := false; two[i]`TestInWhichProved := "hand-set";
end for;
assert not NonHyperellipticAtPrimeCertificate(two[1], 5);
assert NonHyperellipticAtPrimeCertificate(two[2], 3); // #X_0^*(430)(F_9) = 24 > 2*9+2 = 20
assert not NonHyperellipticAtPrimeCertificate(two[3], 7);
assert [X`N : X in two[1..2]] eq [258, 430]; // 258 must precede 430, or it is never recorded and the filter goes untested
uncertified := [];
HHProposition1(~two, ~uncertified);
assert assigned two[4]`IsSubhyp and not two[4]`IsSubhyp and two[4]`TestInWhichProved eq "HHproposition1 isomorphic to 2";
assert not assigned two[5]`IsSubhyp;
assert uncertified eq [<3, 5, 7>];
