// tests/CRV_15_4.m
//
// THE FULL-CURVE POINT-COUNT CHECK THAT `ModelChecks` CANNOT DO.
//
// ⚠ A REAL GAP IN OUR VALIDATION, found 2026-09-06. `VerifyModelSet` skips every `CRV` entry --
//     if Type(e[2]) eq MonStgElt then continue; end if;   // "CRV": non-hyperelliptic entry
// -- so for the paired presentations (`15_4`, `93_1`, `57_1`, `21_2`) the `W={1}` FULL CURVE has
// never been checked against anything except its own recorded genus. Every one of those bases is
// validated only through its hyperelliptic quotients.
//
// That gap had a concrete cost: `X_0^15(4)`'s conic constant `b` is INVISIBLE to the quotient
// checks. `VerifyModelSet` returns 16 checks / 0 failures for EVERY value of `b`, because the
// conic is genus 0 and point counts do not constrain its twist. A test that passes for every
// value of the thing it is supposed to pin is vacuous -- the repo's own recurring lesson.
//
// WHAT THIS TEST DOES. The full curve is the fibre product
//     y^2 = f(x) = -(4x^2-x+1)(4x^2+x+1)(5x^2+3),      z^2 = b*(3x^2+1),
// so #X(F_p) = sum over x in P^1(F_p) of (#y)(#z). Compare that with `ComputePointsViaTrace`,
// which evaluates the Eichler-Selberg trace formula on the W-fixed part of the D-new space and is
// INDEPENDENT of the Borcherds/Schofer machinery -- and, for `15_4`, independent of Guo-Yang too,
// which matters because that base is literature-derived (see data/models/models_15_4.m).
//
// ⚠ IT MUST DISCRIMINATE, OR IT PROVES NOTHING. The test therefore also checks that the WRONG
// constants FAIL. Measured: b = -1 matches at all 12 primes 7..47, while
// b = -2,-3,-5,-6,-7,-10,-15,-30,1,2,3,5,15 each fail at 4-6 of them. This is what closed the
// open `b` question left by 59de486.
//
// ⚠ SCOPE: `15_4` only. Generalising to the other CRV bases means parsing their stored equation
// strings, whose variable names differ per file (`93_1` uses s,z; `15_4` uses x,s), so a general
// version needs the ambient weights recorded in the model files -- the same gap `PROVENANCE.md`
// notes for `ModelRegen` and `21_2`/`57_1`. Worth doing; not done here.

printf "CRV_15_4.m: full-curve point counts vs the trace formula...";

crv_D := 15; crv_N := 4;
crv_X1 := [c : c in GetHyperellipticCandidates()
             | c`D eq crv_D and c`N eq crv_N and c`W eq {1}];
error if #crv_X1 ne 1,
    Sprintf("CRV_15_4: expected exactly one W={1} candidate, found %o", #crv_X1);
crv_X1 := crv_X1[1];

crv_P<crv_x> := PolynomialRing(Rationals());
crv_f := -(4*crv_x^2 - crv_x + 1)*(4*crv_x^2 + crv_x + 1)*(5*crv_x^2 + 3);
crv_g := 3*crv_x^2 + 1;

// #X(F_p) for the fibre product y^2 = f, z^2 = b*g, including the fibre over x = oo.
// f has even degree 6 (y weight 3) and g degree 2 (z weight 1), so both are unramified at oo and
// the count there is governed by the leading coefficients.
function crv_count(b, p, f, g)
    Fp := GF(p);
    Px := PolynomialRing(Fp);
    fp := Px!f; gp := Px!g;
    n := 0;
    for a in Fp do
        fa := Evaluate(fp, a); ga := Fp!b*Evaluate(gp, a);
        ny := (fa eq 0) select 1 else (IsSquare(fa) select 2 else 0);
        nz := (ga eq 0) select 1 else (IsSquare(ga) select 2 else 0);
        n +:= ny*nz;
    end for;
    ny := IsSquare(LeadingCoefficient(fp)) select 2 else 0;
    nz := IsSquare(Fp!b*LeadingCoefficient(gp)) select 2 else 0;
    return n + ny*nz;
end function;

crv_ps := [p : p in [7,11,13,17,19,23,29,31,37,41,43,47] | (crv_D*crv_N) mod p ne 0];
crv_exp := [ComputePointsViaTrace(crv_X1, p, 1) : p in crv_ps];

// (1) the committed constant must MATCH at every prime
crv_nchk := 0; crv_bad := [];
for i->p in crv_ps do
    crv_nchk +:= 1;
    n := crv_count(-1, p, crv_f, crv_g);
    if n ne crv_exp[i] then
        Append(~crv_bad, Sprintf("p=%o: counted %o, trace formula says %o", p, n, crv_exp[i]));
    end if;
end for;
error if not IsEmpty(crv_bad),
    Sprintf("CRV_15_4: the committed b = -1 DISAGREES with the trace formula: %o", crv_bad);

// (2) NEGATIVE CONTROL -- the wrong constants must FAIL, or (1) is vacuous.
crv_wrong := [-2,-3,-5,-6,-7,-10,-15,-30,1,2,3,5,15];
crv_undetected := [];
for b in crv_wrong do
    fails := [p : i->p in crv_ps | crv_count(b, p, crv_f, crv_g) ne crv_exp[i]];
    crv_nchk +:= 1;
    if IsEmpty(fails) then Append(~crv_undetected, b); end if;
end for;
error if not IsEmpty(crv_undetected),
    Sprintf("CRV_15_4: these WRONG conic constants were NOT detected, so the test is vacuous "
            * "for b and models_15_4.m's claim that b is pinned must be weakened: %o",
            crv_undetected);

printf " ok (%o checks over %o primes; b = -1 matches, all %o alternatives rejected)\n",
       crv_nchk, #crv_ps, #crv_wrong;
