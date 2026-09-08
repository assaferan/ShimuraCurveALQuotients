// tests/_rebaselever.m -- NOT a test (leading underscore: excluded from the CI matrix). Run it by
// hand on any base with EMPTY cover keys:
//     magma -b Dd:=10 Nn:=19 tests/_rebaselever.m < /dev/null > log 2>&1
//
// ⚠ AN EMPTY COVER KEY IS NOT EVIDENCE OF AN OBSTRUCTION. The pipeline never builds a genus-g
// curve directly; it assembles it as a FIBRE PRODUCT and needs another quotient to carry an
// equation of degree exactly g+1 over a base shared with the P1/conic. Which degree a quotient's
// equation has depends on whether infinity is a branch point -- that is the HAUPTMODUL
// NORMALISATION, ours to choose, not a fact about the curve. So sweep it before concluding
// anything: t -> r + 1/u at a RATIONAL ROOT r moves a branch point to infinity, turning a quartic
// into a cubic and collapsing a conic that shares the root to degree 1 (a P1).
//
// ⚠⚠ FILLING A KEY IS NOT BUILDING THE RIGHT CURVE, AND THIS TOOL DOES BOTH. Read this before
// trusting any output.
//   * At 22_5 it is CORRECT: it builds the full genus-5 curve and reproduces Guo-Yang's degree-12
//     polynomial VERBATIM, coefficient for coefficient (tests/FullCurve_22_5.m).
//   * At 10_19 it is WRONG. It fills 3 of the 4 empty keys at r = 0 and a different 3 at r = 32/27,
//     but the genus-5 pair it produces has a y-side that is NOT ISOMORPHIC TO ANY Atkin-Lehner
//     quotient of X_0^10(19) -- checked against all three genus-2 quotients derived from Guo-Yang's
//     own curve and involutions (tests/GuoYangQuotients_10_19.m). Its conic IS right, which is
//     exactly why the output looks plausible.
// ⇒ SO NEVER ACCEPT A FILLED KEY WITHOUT AN INDEPENDENT ORACLE. Why it succeeds at 22_5 and fails
// at 10_19 is NOT diagnosed; the likely suspect is that the rebase changes the coordinate but the
// y^2-scale (twist) is not re-derived from CM values afterwards, so the twist can come out wrong.
// ⚠ The root matters and no single one wins -- sweep them all.
//
// ⚠ SLOW on big bases: the AllEquationsAboveCovers call dominates (~40 s at 22_5, ~35 min at
// 10_19). Magma buffers stdout to a file, so an empty log is NOT evidence of a stall -- check CPU
// time on the magma.exe pid, not on the shell wrapper (a `pgrep | head -1` picks the wrapper,
// which sits at 0% and looks exactly like a hang).
//
// GENERAL rebase lever: sweep rational roots on the star base, re-propagate (P1, conic AND
// pointless-conic), and report which previously-EMPTY cover keys get filled.
AttachSpec("ShimuraQuotients.spec");
_ := ClassNumberLU(-4);
import "EquationsCovers.m" : equation_above_P1, equation_above_conic;

D := StringToInteger(Dd); N := StringToInteger(Nn);
curves := GetHyperellipticCandidates();
assert exists(Xstar){X : X in curves | X`D eq D and X`N eq N and IsStarCurve(X)};
covers, ws := AllEquationsAboveCovers(Xstar, curves);

bases := {}; for k in Keys(covers) do bases join:= Keys(covers[k]); end for;
assert exists(STAR){b : b in bases | #{k : k in Keys(covers) | IsDefined(covers[k], b)} ge 3};
printf "star base %o\n", STAR;

baseline := {k : k in Keys(covers) | not IsEmpty(Keys(covers[k]))};
empty0   := {k : k in Keys(covers) | IsEmpty(Keys(covers[k]))};
printf "baseline: %o populated, %o EMPTY: %o\n", #baseline, #empty0,
       Sort([Sprint(Sort(SetToSequence(curves[k]`W))) : k in empty0]);

function rebase_poly(f, r)
    d := Degree(f); e := 2*Ceiling(d/2);
    R<u> := PolynomialRing(Rationals());
    num := &+[ Coefficient(f,i) * (r*u + 1)^i * u^(e-i) : i in [0..d] ];
    while (Degree(num) ge 2) and (Coefficient(num,0) eq 0) and (Coefficient(num,1) eq 0) do
        num := num div u^2; end while;
    return num;
end function;

roots := {Rationals()|};
for k in Keys(covers) do
    if not IsDefined(covers[k], STAR) then continue; end if;
    C := covers[k][STAR];
    if Type(C) ne CrvHyp then continue; end if;
    for rt in Roots(HyperellipticPolynomials(C)) do Include(~roots, rt[1]); end for;
end for;
printf "candidate roots on the star base: %o\n", Sort(Setseq(roots));

best := {}; bestr := 0;
for r in Sort(Setseq(roots)) do
    eqs := AssociativeArray();
    for k in Keys(covers) do
        if not IsDefined(covers[k], STAR) then continue; end if;
        C := covers[k][STAR];
        if Type(C) ne CrvHyp then continue; end if;
        eqs[k] := AssociativeArray();
        eqs[k][STAR] := HyperellipticCurve(rebase_poly(HyperellipticPolynomials(C), r));
    end for;
    changed := true;
    while changed do
        changed := false;
        for L in [1..#curves] do
            if IsDefined(eqs, L) and not IsEmpty(Keys(eqs[L])) then continue; end if;
            if not assigned curves[L]`Covers then continue; end if;
            g := curves[L]`g;
            for Pl in curves[L]`Covers do
                if not IsDefined(eqs, Pl) then continue; end if;
                for b in Keys(eqs[Pl]) do
                    if Type(eqs[Pl][b]) ne CrvHyp then continue; end if;
                    fP := HyperellipticPolynomials(eqs[Pl][b]);
                    isP1  := Degree(fP) eq 1;
                    isCon := (Degree(fP) eq 2) and HasRationalPoint(Conic(eqs[Pl][b]));
                    isPtl := (Degree(fP) eq 2) and not HasRationalPoint(Conic(eqs[Pl][b]));
                    if not (isP1 or isCon or isPtl) then continue; end if;
                    for Ol in curves[L]`Covers do
                        if Ol eq Pl or not IsDefined(eqs, Ol) then continue; end if;
                        if not IsDefined(eqs[Ol], b) then continue; end if;
                        if Type(eqs[Ol][b]) ne CrvHyp then continue; end if;
                        if Degree(HyperellipticPolynomials(eqs[Ol][b])) ne g+1 then continue; end if;
                        if isPtl then
                            wt := curves[Ol]`g + 1;
                            P3<xx,yy,ss,zz> := WeightedProjectiveSpace(Rationals(),[1,wt,1,1]);
                            e1 := Homogenization(Evaluate(fP, ss), zz);
                            e2 := Homogenization(Evaluate(HyperellipticPolynomials(eqs[Ol][b]), ss), zz);
                            H := Curve(P3, [yy^2 - e2, xx^2 - e1]);
                        else
                            H := isP1 select equation_above_P1(eqs[Ol][b], eqs[Pl][b])
                                        else equation_above_conic(eqs[Ol][b], eqs[Pl][b]);
                        end if;
                        if not IsDefined(eqs, L) then eqs[L] := AssociativeArray(); end if;
                        eqs[L][Pl] := H; changed := true; break;
                    end for;
                    if IsDefined(eqs, L) and not IsEmpty(Keys(eqs[L])) then break; end if;
                end for;
                if IsDefined(eqs, L) and not IsEmpty(Keys(eqs[L])) then break; end if;
            end for;
        end for;
    end while;
    got := {k : k in empty0 | IsDefined(eqs, k) and not IsEmpty(Keys(eqs[k]))};
    printf "  r = %-8o fills %o of the %o empty keys: %o\n", r, #got, #empty0,
           Sort([Sprint(Sort(SetToSequence(curves[k]`W))) : k in got]);
    if #got gt #best then best := got; bestr := r; end if;
end for;
printf "\nBEST r = %o filling %o of %o previously-empty keys\n", bestr, #best, #empty0;
exit;
