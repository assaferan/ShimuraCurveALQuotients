// The star of X_0^15(38): a_p from the trace formula for every good prime up to the Sturm bound of level 570
// (240), against 570e4; and the Atkin-Lehner signs against the bad-prime a_p of 570e4 (JL: the sign at p | D flips).
AttachSpec("ShimuraQuotients.spec");
D := 15; N := 38;
W := {m : m in Divisors(D*N) | GCD(m, (D*N) div m) eq 1};
Xstar := CreateShimuraQuot(D, N, W); Xstar`g := GenusShimuraCurveQuotient(D, N, W); Xstar`CurveID := 0;
printf "genus %o\n", Xstar`g;
E := EllipticCurve("570e4");
bad := [];
for p in PrimesInInterval(2, 240) do
    if 570 mod p eq 0 then continue; end if;
    ap := p + 1 - ComputePointsViaTrace(Xstar, p, 1);
    if ap ne TraceOfFrobenius(E, p) then Append(~bad, <p, ap, TraceOfFrobenius(E, p)>); end if;
end for;
printf "good primes <= 240: mismatches with 570e4: %o\n", bad;
printf "570e4 a_p at the bad primes 2, 3, 5, 19: %o\n", [TraceOfFrobenius(E, p) : p in [2,3,5,19]];
curves := GetQuotientsAndGenera([Xstar]); _ := exists(star){c : c in curves | IsStarCurve(c)};
printf "Sturm bound for weight 2, level 570: %o\n", Floor(2 * 570 * (3/2) * (4/3) * (6/5) * (20/19) / 12);
quit;
