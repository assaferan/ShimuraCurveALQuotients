// tests/_offline/FullCurve_22_5.m -- OFFLINE (~420 s; `run_tests.m` globs only `tests/*.m`, so
// `_offline` is excluded automatically). The CI-side check for this base is
// tests/GuoYangQuotients_22_5.m, which takes 0.06 s and validates the STORED entries against
// Guo-Yang; what lives here is the CONSTRUCTION, which needs a full pipeline run.
//
// constructs the FULL curve X_0^22(5), which the pipeline's standard
// stages do not produce, and checks it against Guo-Yang's published equation.
//
// ⚠ WHY THIS EXISTS. models_22_5.m has EMPTY entries at W = {1}, {1,2}, {1,5}, {1,11}: the
// pipeline never builds a genus-g curve directly, it assembles it as a FIBRE PRODUCT, and
// process_P1_cover needs some other quotient to carry an equation of degree exactly g+1 over a
// base shared with the P1/conic. At 22_5 that fails twice and the failures cascade -- the genus-2
// quotients need degree 3 and no degree-3 equation exists (degrees produced: 1,2,4,6,7,8), and
// W={1} needs degree 6 over a base shared with the P1 [1,110], while the only degree-6 equation
// lives over the star base. Guo-Yang have no such trouble because they do not use a fibre product.
//
// ⚠ THE FIX IS A CHANGE OF HAUPTMODUL ON THE STAR BASE, and it is cheap. Which degree a quotient's
// equation has depends on whether the point at infinity is a branch point, i.e. on the Hauptmodul
// normalisation. The Mobius change t -> r + 1/u, for r a RATIONAL ROOT of an equation over that
// base, moves a branch point to infinity: a quartic becomes a cubic, and a conic sharing that root
// collapses to degree 1, a P1. At r = 4/5 both happen at once for {1,2,11,22} (4 -> 3) and
// {1,10,11,110} (2 -> 1), which are exactly the two covers {1,11} needs -- and once the genus-2
// quotients exist, W={1} follows.
//
// ⚠ WHAT MAKES THIS EVIDENCE AND NOT A COINCIDENCE: the result is Guo-Yang's polynomial VERBATIM,
// coefficient for coefficient, not merely isomorphic to it. Nothing in the construction knows
// their equation. r = 0 gives a different model that is isomorphic but not equal, so the test
// pins both: equality at r = 4/5 AND isomorphism at r = 0.

_ := ClassNumberLU(-4);                    // force package loads before `import` (CLAUDE.md)
import "EquationsCovers.m" : equation_above_P1, equation_above_conic;

// y^2 = f(t) with t = r + 1/u becomes Y^2 = u^e f(r + 1/u), Y = y*u^(e/2), e = 2*ceil(deg f/2).
function fc225_rebase(f, r)
    d := Degree(f); e := 2*Ceiling(d/2);
    R<u> := PolynomialRing(Rationals());
    num := &+[ Coefficient(f,i) * (r*u + 1)^i * u^(e-i) : i in [0..d] ];
    while (Degree(num) ge 2) and (Coefficient(num,0) eq 0) and (Coefficient(num,1) eq 0) do
        num := num div u^2;
    end while;
    return num;
end function;

fc225_curves := GetHyperellipticCandidates();
assert exists(fc225_Xstar){X : X in fc225_curves | X`D eq 22 and X`N eq 5 and IsStarCurve(X)};
fc225_covers, fc225_ws := AllEquationsAboveCovers(fc225_Xstar, fc225_curves);

fc225_P<xx> := PolynomialRing(Rationals());
fc225_gy := -11*xx^12 - 80*xx^10 - 240*xx^8 - 362*xx^6 - 240*xx^4 - 80*xx^2 - 11;

// the base every first-level equation sits over
fc225_bases := {};
for k in Keys(fc225_covers) do fc225_bases join:= Keys(fc225_covers[k]); end for;
assert exists(fc225_STAR){b : b in fc225_bases |
        #{k : k in Keys(fc225_covers) | IsDefined(fc225_covers[k], b)} ge 5};

function fc225_build(covers, curves, STAR, r)
    eqs := AssociativeArray();
    for k in Keys(covers) do
        if not IsDefined(covers[k], STAR) then continue; end if;
        C := covers[k][STAR];
        if Type(C) ne CrvHyp then continue; end if;
        eqs[k] := AssociativeArray();
        eqs[k][STAR] := HyperellipticCurve(fc225_rebase(HyperellipticPolynomials(C), r));
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
                    fP := HyperellipticPolynomials(eqs[Pl][b]);
                    isP1  := Degree(fP) eq 1;
                    isCon := (Degree(fP) eq 2) and HasRationalPoint(Conic(eqs[Pl][b]));
                    if not (isP1 or isCon) then continue; end if;
                    for Ol in curves[L]`Covers do
                        if Ol eq Pl or not IsDefined(eqs, Ol) then continue; end if;
                        if not IsDefined(eqs[Ol], b) then continue; end if;
                        if Degree(HyperellipticPolynomials(eqs[Ol][b])) ne g+1 then continue; end if;
                        H := isP1 select equation_above_P1(eqs[Ol][b], eqs[Pl][b])
                                    else equation_above_conic(eqs[Ol][b], eqs[Pl][b]);
                        if not IsDefined(eqs, L) then eqs[L] := AssociativeArray(); end if;
                        eqs[L][Pl] := H; changed := true;
                        break;
                    end for;
                    if IsDefined(eqs, L) and not IsEmpty(Keys(eqs[L])) then break; end if;
                end for;
                if IsDefined(eqs, L) and not IsEmpty(Keys(eqs[L])) then break; end if;
            end for;
        end for;
    end while;
    return eqs;
end function;

assert exists(fc225_l1){k : k in [1..#fc225_curves] | assigned fc225_curves[k]`W
        and fc225_curves[k]`W eq {1} and fc225_curves[k]`D eq 22 and fc225_curves[k]`N eq 5};

fc225_n := 0;

// Each build costs ~90 s, so compute the two roots ONCE and reuse them across (a), (b), (c).
fc225_cache := AssociativeArray();
for fc225_r in [Rationals()| 4/5, 0] do
    fc225_cache[fc225_r] := fc225_build(fc225_covers, fc225_curves, fc225_STAR, fc225_r);
end for;

// (a) r = 4/5 must reproduce Guo-Yang's polynomial EXACTLY.
fc225_e := fc225_cache[4/5];
error if not IsDefined(fc225_e, fc225_l1) or IsEmpty(Keys(fc225_e[fc225_l1])),
    "X0^22(5): the base change at r = 4/5 did not build the full curve W={1}";
for fc225_b in Keys(fc225_e[fc225_l1]) do
    fc225_C := fc225_e[fc225_l1][fc225_b];
    error if Genus(fc225_C) ne 5,
        Sprintf("X0^22(5): built W={1} has genus %o, expected 5", Genus(fc225_C));
    error if HyperellipticPolynomials(fc225_C) ne fc225_gy,
        Sprintf("X0^22(5): built W={1} is NOT Guo-Yang's polynomial verbatim.\n  got %o\n  want %o",
                HyperellipticPolynomials(fc225_C), fc225_gy);
    fc225_n +:= 1;
end for;

// (b) the three genus-2 quotients that were empty must now exist, at r = 4/5 or r = 0.
fc225_g2 := 0;
for fc225_r in [Rationals()| 4/5, 0] do
    fc225_ee := fc225_cache[fc225_r];
    for fc225_W in [{1,2},{1,5},{1,11}] do
        assert exists(fc225_L){k : k in [1..#fc225_curves] | assigned fc225_curves[k]`W
                and fc225_curves[k]`W eq fc225_W and fc225_curves[k]`D eq 22 and fc225_curves[k]`N eq 5};
        if IsDefined(fc225_ee, fc225_L) and not IsEmpty(Keys(fc225_ee[fc225_L])) then
            for fc225_b in Keys(fc225_ee[fc225_L]) do
                error if Genus(fc225_ee[fc225_L][fc225_b]) ne 2,
                    Sprintf("X0^22(5): built W=%o has genus %o, expected 2",
                            fc225_W, Genus(fc225_ee[fc225_L][fc225_b]));
            end for;
            fc225_g2 +:= 1;
        end if;
    end for;
end for;
error if fc225_g2 lt 3,
    Sprintf("X0^22(5): expected all three genus-2 quotients to be built, got %o", fc225_g2);

// (c) r = 0 must give an ISOMORPHIC but NOT equal model -- so (a) is pinning a real coincidence
// of coordinates, not just any construction that happens to land on genus 5.
fc225_e0 := fc225_cache[0];
error if not IsDefined(fc225_e0, fc225_l1) or IsEmpty(Keys(fc225_e0[fc225_l1])),
    "X0^22(5): the base change at r = 0 did not build the full curve";
for fc225_b in Keys(fc225_e0[fc225_l1]) do
    fc225_C0 := fc225_e0[fc225_l1][fc225_b];
    error if not IsIsomorphic(fc225_C0, HyperellipticCurve(fc225_gy)),
        "X0^22(5): the r = 0 model is not isomorphic to Guo-Yang's";
    error if HyperellipticPolynomials(fc225_C0) eq fc225_gy,
        "X0^22(5): r = 0 gave Guo-Yang's polynomial verbatim too -- the r = 4/5 equality is then "
        * "not evidence of anything, and this test's claim must be weakened";
    fc225_n +:= 1;
end for;

printf " ok (X0^22(5) full curve: %o comparison(s); Guo-Yang's degree-12 polynomial reproduced "
       * "VERBATIM at r = 4/5, %o genus-2 quotients recovered)\n", fc225_n, fc225_g2;
