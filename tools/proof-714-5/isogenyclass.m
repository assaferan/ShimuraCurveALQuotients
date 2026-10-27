// a_p of a star against a Cremona curve for every good prime up to the Sturm bound of its level.
//   magma -b DD:=6 NN:=385 LAB:=2310q1 sturm.m
AttachSpec("ShimuraQuotients.spec");
D := StringToInteger(DD); N := StringToInteger(NN);
W := {m : m in Divisors(D*N) | GCD(m, (D*N) div m) eq 1};
Xstar := CreateShimuraQuot(D, N, W); Xstar`g := GenusShimuraCurveQuotient(D, N, W); Xstar`CurveID := 0;
E := EllipticCurve(LAB); lev := Conductor(E);
bound := Floor(2 * lev * &*[Rationals() | 1 + 1/q : q in PrimeDivisors(lev)] / 12);
printf "%o_%o vs %o: Sturm bound %o\n", D, N, LAB, bound;
t := Realtime(); bad := [];
for p in [q : q in PrimesInInterval(2, bound) | lev mod q ne 0] do
    ap := p + 1 - ComputePointsViaTrace(Xstar, p, 1);
    if ap ne TraceOfFrobenius(E, p) then Append(~bad, p); end if;
end for;
printf "mismatches: %o (%o s)\n", bad, Realtime() - t;
quit;
