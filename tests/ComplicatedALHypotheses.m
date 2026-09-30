// tests/ComplicatedALHypotheses.m
//
// CheckComplicatedALHypotheses checks the hypotheses of prop:complicatedAL for
// TestComplicatedALFixedPointsOnQuotient (G an AL group) and CheckGeneralizedComplicatedFixedPoints
// (G = <W_odd, V2>).  This pins:
//   * [FH99] Table 3 (p. 116): the 34 quotients X_0(N)/W' that [FH] Prop 6 proves non-hyperelliptic
//     are still proved (VerifyFHTable3), with the witness of X_0(58)/<w2>;
//   * one negative control for each rejection reason of CheckComplicatedALHypotheses;
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
// X_0(58)/<w2>: h(-116) = 6 (LMFDB 2.0.116.1); genus(X_0(58)/<w2>) = 3 (Magma modular symbols).
// The third entry of the witness depends on the iteration order of W, so it is not pinned.
cah_S58 := CuspidalSubspace(ModularSymbols(58, 2, 1));
assert Dimension(Kernel(AtkinLehner(cah_S58, 2) - 1)) eq 3 and ClassNumber(-116) eq 6;
cah_58 := TestComplicatedALFixedPointsOnQuotient(1, 58);
assert cah_58[{1, 2}][1..2] eq [58, 29];
assert CheckComplicatedALHypotheses(1, 58, {1, 2}, 58, 29);

// One negative control per rejection reason; each case fails only the named hypothesis among
// those checked before it.
cah_reason := func<D, N, W, N1, N2, V | why where _, why := CheckComplicatedALHypotheses(D, N, W,
                                                                               N1, N2 : V := V)>;
assert cah_reason(1, 58, {1, 2}, 29, 29, 0) eq "N1 N2";          // N1 = N2
assert cah_reason(1, 28, {1, 4}, 2, 7, 0) eq "N1 N2";            // 2 | 28/2: w2 is not an AL involution
assert cah_reason(1, 58, {1, 2}, 3, 29, 0) eq "divisor";         // 3 does not divide 58
assert cah_reason(1, 28, {1, 4}, 4, 7, 0) eq "congruence";       // 7 = 7 mod 8 and D = 1 is odd
// h(-68) = 4 (LMFDB 2.0.68.1).
assert cah_reason(1, 136, {1, 136}, 8, 17, 0) eq "class number";
// nu(w2) on X_0(22) is 2, not #G = 1: genus(X_0(22)) = 2, genus(X_0(22)/<w2>) = 1 (Magma
// modular symbols), and Riemann-Hurwitz 2 = 2 (2 - 2) + nu.
cah_S22 := CuspidalSubspace(ModularSymbols(22, 2, 1));
assert Dimension(cah_S22) eq 2 and Dimension(Kernel(AtkinLehner(cah_S22, 2) - 1)) eq 1;
assert cah_reason(1, 22, {1}, 2, 11, 0) eq "nu(N1)";
assert cah_reason(1, 58, {2, 58}, 58, 29, 0) eq "group";         // 1 is not in W
assert cah_reason(1, 58, {1, 29}, 58, 29, 0) eq "notin G";       // 29 is in W
// V3 on X_0(198) conjugates w_m to w_{9m} for m = 2 mod 3 (see CheckGeneralizedComplicatedFixedPoints).
cah_V3 := get_V3(198);
assert cah_reason(1, 198, {1, 2}, 2, 11, cah_V3) eq "group";     // V3 and w2 do not commute
assert cah_reason(1, 198, {1, 22}, 2, 99, cah_V3) eq "product in G";  // w2 w99 = w198, not in {w1, w22}
assert cah_reason(1, 198, {1, 22}, 2, 11, cah_V3) eq "commute";  // V3 w11 V3^-1 = w99

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
