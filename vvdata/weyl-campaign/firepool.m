// How many FIRING DEGREE-1 discriminants does a base actually have?
//
// A row of the HauptmodulM0Residuals system is INFORMATIVE only if its firing status
// differs from a normaliser's, and rows with degs ne 1 are discarded outright.  So the
// sweep needs >= 2 firing degree-1 discriminants beyond the normalisers to pin (r1,r2).
// 22_3 turned out to have only TWO in total, one of them already an anchor -- if the stuck
// bases (33_2, 34_3, 46_3) are similarly starved then the FIRE knob has nothing to give and
// they are undeterminable by this sweep for want of CM points on the base itself.
//
// bd is the class-number bound RationalandQuadraticCMPoints scans to (default 4).  Raising
// it is the other lever, so report both.
//
//   magma -b firepool.m
AttachSpec("ShimuraQuotients.spec");

bases := [[33,2],[34,3],[46,3],[22,3],[38,3],[21,2],[15,2]];
BD := 4; if assigned BDD then BD := StringToInteger(BDD); end if;

for b in bases do
    D := b[1]; N := b[2];
    Xstar := CreateShimuraQuot(D, N, Set(Divisors(D*N)));
    Xstar`g := GenusShimuraCurveQuotient(D, N, Xstar`W); Xstar`CurveID := 0;
    curves := GetQuotientsAndGenera([Xstar]);
    _ := exists(star){c : c in curves | IsStarCurve(c)};
    rat_all := RationalandQuadraticCMPoints(star : coprime_to_level := false, bd := BD);
    fires := func<d | not IsEmpty(PrimeDivisors(N div GCD(N, FundamentalDiscriminant(d))))>;
    f1 := [ p[1] : p in rat_all | fires(p[1]) ];
    n1 := [ p[1] : p in rat_all | not fires(p[1]) ];
    printf "X0^%o(%o) bd=%o: %o rational pts; FIRING deg-1 %o ; non-firing deg-1 %o\n",
           D, N, BD, #rat_all, f1, n1;
end for;
printf "DONE\n";
quit;
