// tests/GeneralizedComplicatedV3.m
//
// FilterByGeneralizedComplicatedFixedPoints (the generalized [FH] Prop 6, prop:complicatedAL)
// needs G = <W_odd, V_p> to be an elementary 2-group with X/W = X/G and w_N1, w_N2 notin G.
//
// V3 branch (removed).  The old commutation check (bad_sets) only ran when 9 notin W, but the
// certificate requires 9 in W, so it never ran.  All four recorded V3 certificates had
// V3 w_m V3^-1 = w_{9m} for m = 2 mod 3, so G contained w_9, was non-abelian of order 16, and
// N1, N2 lay in G: the certificates were invalid.  When every w in W_odd commutes with V3,
// Fix(V3 on C) >= 8 contradicts the guard Fix = 4, so no V3 certificate can ever be valid.
//
// This pins:
//   * the four former V3 certificates (CurveIDs 5124, 7923, 8387, 9255) no longer fire;
//   * the commuting V3 control (5,63)/<7,9> does not fire;
//   * the valid V2 certificate #997, X_0(216)/<w8,w27> (G = <w27, V2>, N1 = 8, N2 = 216),
//     still fires, so the new guard (S2 commutes with W_odd, G elementary abelian of order #W,
//     w_N1, w_N2 notin G) does not reject a genuine certificate.

gcv3_make := function(D, N, gens, id)
    X := CreateShimuraQuot(D, N, AllALsFromGens(gens, D*N));
    X`g := GenusShimuraCurveQuotient(D, N, X`W);
    X`CurveID := id;
    return X;
end function;

gcv3_bad := [<10, 153, {2, 9, 85}, 5124>, <22, 45, {2, 9, 55}, 7923>,
             <26, 45, {5, 9, 26}, 8387>, <35, 18, {2, 9, 35}, 9255>,
             <5, 63, {7, 9}, 0>];
gcv3_curves := [gcv3_make(t[1], t[2], t[3], t[4]) : t in gcv3_bad];
for gcv3_X in gcv3_curves do
    gcv3_ok := CheckGeneralizedComplicatedFixedPoints(gcv3_X);
    assert not gcv3_ok;
end for;
FilterByGeneralizedComplicatedFixedPoints(~gcv3_curves);
assert &and[not assigned c`IsSubhyp and not assigned c`IsHyp : c in gcv3_curves];

gcv3_997 := gcv3_make(1, 216, {8, 27}, 997);
assert gcv3_997`g eq 5;
gcv3_ok, gcv3_s := CheckGeneralizedComplicatedFixedPoints(gcv3_997);
assert gcv3_ok;
assert gcv3_s eq "GeneralizedComplicatedFixedPoints: V2, pPart=8, N1=8 (nu=4), N2=216 (nu=12)";
gcv3_list := [gcv3_997];
FilterByGeneralizedComplicatedFixedPoints(~gcv3_list);
assert gcv3_list[1]`IsHyp eq false and gcv3_list[1]`IsSubhyp eq false;
assert gcv3_list[1]`TestInWhichProved eq gcv3_s;
