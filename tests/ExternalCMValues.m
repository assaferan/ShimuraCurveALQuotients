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
// COST, and ADDING MORE CASES.  One base costs a full Borcherds + Schofer run (~850 s here), comparable
// to the X0_D_N model tests.  More published tables are available in the same paper and are easy to add
// to `cases` -- Example 36 gives s(tau_-280) = 5/16 on X_0^14(5) (another odd N, with 5 | 280), and
// Examples 35, 42 and 43 tabulate values on X_0^146(1), X_0^142(1) and X_0^302(1).  Each added case
// costs another such run, so they are left out deliberately rather than overlooked.

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

        // the point of the test: at least one target must be divisible by N, or the row rescaling
        // absorbs the check and it proves nothing
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
