// External validation of Schofer CM values against PUBLISHED singular moduli.
//
// WHY THIS TEST EXISTS, AND WHY THE EXISTING MODEL TESTS DO NOT COVER IT.
// Every X0_D_N test compares a computed CURVE against a published equation. `ReduceTable` is free to
// rescale each row of the Schofer table, and an equation is insensitive to that rescaling. So an error
// that shifts a whole row -- for instance a wrong coefficient of log p at the level prime, which acts
// uniformly across all the discriminants COPRIME to N -- passes those tests silently. That is not
// hypothetical: it is exactly how the odd-N bases pass on main while a genuine per-cover correction
// goes undetected.
//
// A published value of a Hauptmodul at a CM point pins the ABSOLUTE calibration instead, provided the
// point is chosen so the row scaling cannot absorb it: normalise the Hauptmodul at discriminants that
// are coprime to N, then evaluate at one DIVISIBLE by N. The discriminants divisible by N are exactly
// the ones the pipeline's coprime-to-level filter discards, which is why `Keep` is needed to admit
// them (see the Keep parameter of ValuesAtCMPoints).
//
// COMPARING ACROSS NORMALISATIONS. A Hauptmodul is only defined up to a Mobius transformation, and the
// published normalisation generally differs from the pipeline's. The comparison is therefore: fit the
// unique Mobius map sending the computed values at three reference CM points to the published values
// there, then check the remaining published values are reproduced exactly. Three points fix the map, so
// every further point is a genuine check.
//
// SOURCE. J.-W. Guo and Y. Yang, "Equations of hyperelliptic Shimura curves", Compositio/Proc. LMS
// (doi:10.1112/S0010437X16007739), Section 4.2, Example 37: for X = X_0^10(19), with s the Hauptmodul
// of X/W_{10,19} normalised by s(tau_-8) = 0, s(tau_-40) = infinity, s(tau_-3) = 1, Schofer's formula
// gives s(tau_-760) = 32/5.  Note 19 | 760 while 19 does not divide 8, 40 or 3 -- precisely the
// configuration described above.
//
// SENSITIVITY (checked, not assumed): if the Schofer values at the FIRING discriminants -3, -8, -40
// were multiplied by 19^(+-1) -- the shape of error this test exists to catch -- the reconstructed
// s(tau_-760) becomes -32/481 or 608/581 rather than 32/5.  So the check bites.
//
// COST, and ADDING MORE CASES.  One base costs a full Borcherds + Schofer run (~3300 s as measured in
// the suite), comparable to the X0_D_N model tests.
//
// ⚠⚠ THIS PARAGRAPH USED TO RECOMMEND GUO-YANG EXAMPLE 36 AS AN EASY ADDITION. IT IS NOT AN ADDITION
// AT ALL, and the correction is kept rather than deleted because the recommendation was followed far
// enough to nearly cost an hour of oracle run on 2026-09-24.  A usable case needs all three
// conditions now enforced at the guard below, and Guo-Yang's remaining published values fail them:
//
//     Example 36, X_0^14(5), s(tau_-280) = 5/16   FAILS (b): it normalises at -4, -11, -35 and
//                                                 5 | 35, so the value is absorbed by the Mobius fit
//     Examples 35, 42, 43 on X_0^146(1),          FAIL (a): N = 1, where the coprime-to-level filter
//     X_0^142(1), X_0^302(1)                      discards nothing and no value can escape the row
//                                                 rescaling -- the configuration cannot exist
//
// ⇒ Example 37 (X_0^10(19)) is the ONLY Guo-Yang case that satisfies the configuration, which is why
// this test rests on a single published value.  That is a real limit of the available data, not an
// omission: adding more needs a published CM value on an N > 1 base, normalised at discriminants
// coprime to N.  Do not spend a run before checking (a), (b), (c) -- it costs nothing and it is what
// makes or breaks the test.

// The Mobius map sending z0 -> 0, z1 -> infinity, z2 -> 1, evaluated at z.  Any of the arguments may be
// Infinity.  Returns Infinity when the image is the point at infinity.
function mobius(z0, z1, z2, z)
    // cross-ratio (z, z0; z1, z2) written so each factor can be dropped when its argument is Infinity
    num := (z eq Infinity() or z0 eq Infinity()) select 1 else z - z0;
    den := (z eq Infinity() or z1 eq Infinity()) select 1 else z - z1;
    c1  := (z2 eq Infinity() or z1 eq Infinity()) select 1 else z2 - z1;
    c2  := (z2 eq Infinity() or z0 eq Infinity()) select 1 else z2 - z0;
    if den*c2 eq 0 then return Infinity(); end if;
    return (num*c1)/(den*c2);
end function;

procedure test_ExternalCMValues()
    printf "Testing Schofer CM values against published singular moduli...";

    // <D, N, <d0, d1, d2> normalising to <0, oo, 1>, [<disc, published value>]>
    cases := [* <10, 19, <-8, -40, -3>, [* <-760, Rationals()!(32/5)> *]> *];

    nchecked := 0;
    for cs in cases do
        D, N, norm, targets := Explode(cs);
        d0, d1, d2 := Explode(norm);
        keep := {d0, d1, d2} join {t[1] : t in targets};

        // ⚠⚠ THREE CONDITIONS, AND THE FIRST TWO WERE MISSING -- the guard below used to test only
        // the third, which is VACUOUS AT N = 1 because every integer is divisible by 1. A case added
        // on an N = 1 base would therefore have passed this guard while proving nothing, which is
        // precisely the failure mode this file exists to prevent. Checked 2026-09-24.
        //
        //   (a) N > 1. At N = 1 the coprime-to-level filter discards nothing, so there is no
        //       discriminant whose value escapes the row rescaling, and the configuration this test
        //       depends on cannot exist AT ALL. ⇒ Guo-Yang Examples 35, 42 and 43 (on X_0^146(1),
        //       X_0^142(1), X_0^302(1)) can NEVER serve here, however many published values they
        //       carry.
        //   (b) every NORMALISING discriminant coprime to N. The Mobius map is fitted at those three
        //       points, so if one of them is divisible by N its value is absorbed into the fit and
        //       the calibration is lost. ⇒ Guo-Yang Example 36 (X_0^14(5), s(tau_-280) = 5/16)
        //       CANNOT serve either: it normalises at -4, -11 and -35, and 5 | 35. This file's own
        //       header used to recommend Example 36 as an easy addition. It is not an addition at
        //       all, and an hour of oracle run would have been spent discovering that.
        //       ⚠ Independently refuted: [[gy-example36-cannot-discriminate]] shows that value is
        //       reproduced EXACTLY whether corrected or uncorrected, because -35 and -280 are both
        //       non-firing and any correction cancels in the ratio.
        //   (c) at least one TARGET divisible by N -- the original condition, still necessary.
        error if N le 1,
            Sprintf("X0^%o(%o): N = %o, so 'divisible by N' is vacuous (every integer is) and no "
                    * "value here can escape ReduceTable's per-row rescaling. This case cannot "
                    * "calibrate anything; it needs a base with N > 1.", D, N, N);
        error if exists{d : d in norm | GCD(d, N) ne 1},
            Sprintf("X0^%o(%o): a NORMALISING discriminant is not coprime to N, so its value is "
                    * "absorbed by the Mobius fit and the calibration is lost. Normalise at "
                    * "discriminants coprime to N and evaluate at one divisible by N.", D, N);
        error if not exists{t : t in targets | t[1] mod N eq 0},
            Sprintf("X0^%o(%o): no target discriminant is divisible by N; the check would be vacuous",
                    D, N);

        Xstar := CreateShimuraQuot(D, N, Set(Divisors(D*N)));
        Xstar`g := GenusShimuraCurveQuotient(D, N, Xstar`W); Xstar`CurveID := 0;
        curves := GetQuotientsAndGenera([Xstar]);
        assert exists(star){c : c in curves | IsStarCurve(c)};
        tab := ValuesAtCMPoints(star, curves : Keep := keep);

        discs := tab`Discs;
        srow := tab`Values[tab`sIndex];
        idx := AssociativeArray();
        for i->d in discs do idx[d] := i; end for;
        for d in keep do
            error if not IsDefined(idx, d),
                Sprintf("X0^%o(%o): discriminant %o is missing from the table (Keep did not admit it); "
                        * "table has %o", D, N, d, discs);
        end for;

        z0 := srow[idx[d0]]; z1 := srow[idx[d1]]; z2 := srow[idx[d2]];
        for t in targets do
            d, expected := Explode(t);
            got := mobius(z0, z1, z2, srow[idx[d]]);
            error if got ne expected,
                Sprintf("X0^%o(%o): s(tau_%o) = %o, published value is %o "
                        * "(normalised at %o, %o, %o -> 0, oo, 1)",
                        D, N, d, got, expected, d0, d1, d2);
            nchecked +:= 1;
        end for;
    end for;

    printf " OK (%o published value(s))\n", nchecked;
end procedure;

test_ExternalCMValues();
