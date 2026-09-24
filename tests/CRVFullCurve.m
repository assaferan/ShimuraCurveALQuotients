// tests/CRVFullCurve.m
//
// FULL-CURVE PROOFS FOR THE PAIRED (CRV) PRESENTATIONS, WITHOUT CALLING IsIsomorphic ON THEM.
//
// ⚠ THE PROBLEM. `IsIsomorphic`'s cost depends violently on PRESENTATION, not genus
// (`tests/IsoScreen.m`): a genus-7 HYPERELLIPTIC curve settles in 0.06 s, while the genus-3 CRV
// pair at `14_3` runs >50 min and the genus-5 one at `26_3` >1 h. So the paired bases -- `93_1`,
// `26_3`, `15_4` -- have only ever been checked against Guo-Yang at the level of their QUOTIENTS.
// Their model headers say so explicitly.
//
// ⚠ AND SCREENS ARE NOT A SUBSTITUTE. `ScreenByPlaces` and trace-formula point counts can only
// REFUTE: non-isomorphic curves with isogenous Jacobians agree at every prime. Replacing a proof
// with a screen would weaken the claim, not preserve it.
//
// THE CONSTRUCTION. An isomorphism respecting the labelled involutions descends to the common
// genus-0 base as a MOBIUS map `mu`. So:
//   1. the `y`-quotient `y^2 = f` is HYPERELLIPTIC, where `IsIsomorphic` is fast -- get `mu` there;
//   2. read `mu` off the returned map's defining equations;
//   3. require `f_gy(mu)/f_our` AND `g_gy(mu)/g_our` to be CONSTANT SQUARES, which yields the two
//      scalars;
//   4. build the map explicitly and let `IsIsomorphism` certify it.
// Step 4 is what makes this a PROOF rather than a coefficient coincidence. Measured at `93_1`:
// the whole thing runs in ~0.07 s against hours for the generic call.
//
// ⚠ WHY STEP 5 (the Aut fallback) EXISTS. `IsIsomorphic` returns an ARBITRARY element of
// `Isom` -- a torsor under `Aut` -- so the first `mu` it hands back need not be the one that also
// carries the conic. When it is not, compose with automorphisms of the target quotient and retry.
// Without that, a correct pair can be rejected on an unlucky draw.

printf "Proving CRV full-curve isomorphisms against Guo-Yang (constructed, not searched)...";

// Given the two pairs as polynomials in one variable, return whether an isomorphism was
// CONSTRUCTED, plus mu and the two scalars.
function crv_iso_data(f_our, g_our, f_gy, g_gy)
    Cy_our := HyperellipticCurve(f_our);
    Cy_gy  := HyperellipticCurve(f_gy);
    if Genus(Cy_our) ne Genus(Cy_gy) then return false, 0, 0, 0; end if;
    // ⚠ Magma returns ONE value (false) when the curves are not isomorphic, and TWO when they
    // are -- so `ok, phi := ...` errors outright in the negative case. Test first, then fetch.
    if not IsIsomorphic(Cy_our, Cy_gy) then return false, 0, 0, 0; end if;
    _, phi := IsIsomorphic(Cy_our, Cy_gy);

    Pt<t> := Parent(f_our);
    // candidate maps: phi, then phi composed with the target's automorphisms
    cands := [phi];
    try
        A, mA := AutomorphismGroup(Cy_gy);
        cands := [phi*mA(a) : a in A];
    catch e ; end try;

    for psi in cands do
        de := DefiningEquations(psi);
        if #de lt 3 then continue; end if;
        // the base coordinate transforms as x -> de[1], z -> de[3]; dehomogenise to get mu(t)
        num := Evaluate(de[1], [t, 0, 1]);
        den := Evaluate(de[3], [t, 0, 1]);
        if den eq 0 then continue; end if;
        mu := num/den;
        qf := Evaluate(f_gy, mu) / f_our;
        qg := Evaluate(g_gy, mu) / g_our;
        if not (IsCoercible(Rationals(), qf) and IsCoercible(Rationals(), qg)) then continue; end if;
        sf, rf := IsSquare(Rationals()!qf);
        sg, rg := IsSquare(Rationals()!qg);
        if sf and sg then return true, mu, rf, rg; end if;
    end for;
    return false, 0, 0, 0;
end function;

crv_P<t> := PolynomialRing(Rationals());
// <name, weights, f_our, g_our, f_gy, g_gy>  with y^2 = f and x^2 = g on BOTH sides
crv_cases := [*
    // 93_1: our W={1} against Guo-Yang, reading their `-3t` as `-3s` (see GuoYangEquations.m).
    <"93_1", 3,
     t^6 - 4*t^5 + 50/9*t^4 - 34/9*t^3 + 17/9*t^2 - 2/3*t + 1/9,
     -(144*t^2 - 36*t + 63),
     (3*t^3 - 7*t^2 - 3*t - 1)*(3*t^3 + t^2 - 3*t - 9),
     -4*t^2 - 6*t - 9>,

    // 26_3: the committed W={1} entry is the base_label := 8103 presentation, which is the V_4
    // Guo-Yang use. Our f is theirs scaled by 1/64 = (1/8)^2 and the conic is theirs verbatim, so
    // mu is the identity and the map is (t, y, x) -> (t, 8y, x).
    <"26_3", 3,
     1/64*t^6 - 1/32*t^4 + 9/64*t^2 + 1/8,
     -8*t^2 - 3,
     t^6 - 2*t^4 + 9*t^2 + 8,
     -8*t^2 - 3>
*];

// ⚠ HOW 26_3 GOT HERE, because it did not work at first. Its DEFAULT presentation pairs a
// different V_4 than Guo-Yang's: y-side [1,6] and the third [1,26] conic, versus their [1,78] and
// the second. Measured then: the genus-2 y-curves were NOT isomorphic and the conic discriminants
// differed by a non-square. Both V_4s are legitimate -- IsoScreen.m warns the Klein four-group is
// not unique in Aut(C).
// The fix was assaferan's: a CRV is built over a CHOSEN BASE, and the base decides the V_4. Running
// AllEquationsAboveCovers with base_label := 8103 yields the V_4 Guo-Yang use, after which the
// comparison is nearly the identity. models_26_3.m now stores that presentation.
// ⇒ GENERAL LESSON: when a CRV pair will not match, try another base before concluding anything
// about the curve. The presentation is a choice; the curve is not.

crv_n := 0; crv_fail := [];
for c in crv_cases do
    nm, wy, fo, go, fg, gg := Explode(c);
    got, mu, rf, rg := crv_iso_data(fo, go, fg, gg);
    crv_n +:= 1;
    if not got then
        Append(~crv_fail, nm cat ": no isomorphism could be CONSTRUCTED (mu carrying both sides)");
        continue;
    end if;
    // certify the constructed map on the projective curves -- this is the step that proves it
    Pw<x,y,s,z> := WeightedProjectiveSpace(Rationals(), [1,wy,1,1]);
    hom_ := func<p | &+[Coefficient(p,i)*s^i*z^(Degree(p)-i) : i in [0..Degree(p)]]>;
    C1 := Curve(Pw, [y^2 - hom_(fo), x^2 - hom_(go)]);
    C2 := Curve(Pw, [y^2 - hom_(fg), x^2 - hom_(gg)]);
    munum := Numerator(mu); muden := Denominator(mu);
    msn := hom_(munum); msd := (muden eq 1) select z^Degree(munum) else hom_(muden);
    ok := false;
    try
        // ⚠ rg and rf MULTIPLY: rf^2 = f_gy(mu)/f_our means y -> rf*y, and likewise x -> rg*x.
        // Writing x/rg here inverted the conic scalar and the map failed to certify.
        psi := map< C1 -> C2 | [ rg*x, rf*y, msn, msd ] >;
        ok := IsIsomorphism(psi);
    catch e ok := false; end try;
    if not ok then
        Append(~crv_fail, nm cat ": mu and the scalars were found, but the constructed map did not "
                          * "certify as an isomorphism");
    end if;
end for;

error if crv_n eq 0, "CRVFullCurve: NO EVIDENCE -- zero cases checked.";
error if not IsEmpty(crv_fail),
    Sprintf("CRVFullCurve: %o case(s) failed: %o", #crv_fail, crv_fail);

printf " ok (%o CRV pair(s) proven isomorphic to Guo-Yang by explicit construction)\n", crv_n;
