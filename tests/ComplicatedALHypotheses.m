// tests/ComplicatedALHypotheses.m
//
// CheckComplicatedALHypotheses checks the hypotheses of prop:complicatedAL for
// TestComplicatedALFixedPointsOnQuotient (G an AL group) and CheckGeneralizedComplicatedFixedPoints
// (G = <W_odd, V2>).  This pins:
//   * [FH99] Table 3 (p. 116): the 34 quotients X_0(N)/W' that [FH] Prop 6 proves non-hyperelliptic
//     are still proved (VerifyFHTable3), with the witness of X_0(58)/<w2>;
//   * nu(w_N2) = 9 #G is rejected: X_0(118)/<w118> with N1 = 2, N2 = 59 has nu(w59) = 18 = 9 * 2
//     and meets every other hypothesis;
//   * w_N2 hyperelliptic is rejected: X_0(6,143)/<w6,w22,w26> (genus 3) is hyperelliptic, with
//     hyperelliptic involution w286 (genus(X/<w286>) = 0), yet N1 = 13, N2 = 286 meets every other
//     hypothesis;
//   * the mixed group on X_0^*(216) (G = <w27, V2>, N1 = 8, N2 = 216), with the genus of
//     X_0(216)/<G, w216> by TraceDNewQuotient equal to the Riemann-Hurwitz count;
//   * the hyperelliptic X_0^*(136) (V2 is its hyperelliptic involution) is not proved.

cah_make := function(D, N, W)
    X := CreateShimuraQuot(D, N, W);
    X`g := GenusShimuraCurveQuotient(D, N, W);
    return X;
end function;

// [FH99] Table 3, via the pipeline's own verification on every AL quotient of its levels.
cah_levels := [58, 76, 86, 102, 106, 114, 122, 124, 130, 132, 134, 140, 150, 170, 174, 182, 186,
               190, 198, 204, 210, 222, 230, 330, 390];
cah_curves := [cah_make(1, N, Wt[1]) : Wt in ALSubgroups(N), N in cah_levels];
FilterByComplicatedALFixedPointsOnQuotient(~cah_curves);
VerifyFHTable3(cah_curves);
cah_58 := TestComplicatedALFixedPointsOnQuotient(1, 58);
assert cah_58[{1, 2}] eq [58, 29, 29];
assert CheckComplicatedALHypotheses(1, 58, {1, 2}, 58, 29);

// nu(w_N2) = 9 #G.
assert NumFixedPoints(1, 118, 59) eq 18 and NumFixedPoints(1, 118, 2) eq 2;
cah_ok, cah_why := CheckComplicatedALHypotheses(1, 118, {1, 118}, 2, 59);
assert not cah_ok and cah_why eq "nu(N2)";
assert not IsDefined(TestComplicatedALFixedPointsOnQuotient(1, 118), {1, 118});

// w_N2 is the hyperelliptic involution.
cah_W := AllALsFromGens({6, 22, 26}, 858);
assert GenusShimuraCurveQuotient(6, 143, cah_W) eq 3;
assert GenusShimuraCurveQuotient(6, 143, AllALsFromGens({6, 22, 26, 286}, 858)) eq 0;
cah_ok, cah_why := CheckComplicatedALHypotheses(6, 143, cah_W, 13, 286);
assert not cah_ok and cah_why eq "hyperelliptic N2";
assert not IsDefined(TestComplicatedALFixedPointsOnQuotient(6, 143), cah_W);

// Mixed group on X_0^*(216); the full check is pinned in tests/GeneralizedComplicatedV3.m.
cah_V := get_V2(216);
assert CheckComplicatedALHypotheses(1, 216, {1, 27}, 8, 216 : V := cah_V);
cah_H := {1, 27, 8, 216};
cah_nus := [NumFixedPoints(1, 216, w) : w in cah_H | w ne 1] cat
           [NumFixedPointsNonALOnX(cah_V, "V2", w, 1, 216) : w in cah_H];
cah_gRH := ((2*GenusShimuraCurve(1, 216) - 2 - &+cah_nus)/(2*#cah_H) + 2)/2;
assert TraceDNewQuotient(cah_V, "V2", 1, cah_H, 1, 216) eq 2 and cah_gRH eq 2;

// X_0^*(136) is hyperelliptic.
cah_136 := cah_make(1, 136, {1, 8, 17, 136});
assert cah_136`g eq 3;
assert not CheckGeneralizedComplicatedFixedPoints(cah_136);
assert not IsDefined(TestComplicatedALFixedPointsOnQuotient(1, 136), {1, 8, 17, 136});
