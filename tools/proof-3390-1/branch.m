AttachSpec("ShimuraQuotients.spec");
D := 3390; N := 1;
H := {1, 5, 6, 30, 113, 565, 678, 3390};
G := { m : m in Divisors(D*N) | GCD(m, (D*N) div m) eq 1 };
printf "genus Y %o, star %o\n", GenusShimuraCurveQuotient(D, N, H), GenusShimuraCurveQuotient(D, N, G);
for m in Sort(Setseq(G diff H)) do
    e := NumFixedPointsByCMOrder(D, N, m);
    printf "   w_%o: %o fixed points %o\n", m, NumFixedPoints(D, N, m), [ <d, e[d]> : d in Sort(Setseq(Keys(e))) ];
end for;
Xstar := CreateShimuraQuot(D, N, G); Xstar`g := GenusShimuraCurveQuotient(D, N, G); Xstar`CurveID := 0;
curves := GetQuotientsAndGenera([Xstar]); _ := exists(star){c : c in curves | IsStarCurve(c)};
printf "degrees of CM places on the star: %o\n", [ <d, DegreeOfFieldOfDefinitionOfCMPoint(star, d)> : d in [-3, -40, -43, -67, -123, -168, -132, -312, -148, -372, -280, -555] ];
g := 0; t := Realtime();
for p in [q : q in PrimesInInterval(7, 120) | D mod q ne 0] do
    s1 := ComputePointsViaTrace(Xstar, p, 1); s2 := ComputePointsViaTrace(Xstar, p, 2);
    c1 := s1 - p - 1; c2 := (s2 - p^2 - 1 + c1^2) div 2;
    L1 := 1 + c1 + c2 + p*c1 + p^2;          // #J(F_p) = L_p(1)
    g := GCD(g, L1);
    if p le 13 then printf "p = %o: c1 = %o, c2 = %o, T(p^2) char poly x^2 + %o x + %o, #J(F_p) = %o\n", p, c1, c2, c1, c2 - 2*p, L1; end if;
end for;
printf "torsion bound gcd #J(F_p), p < 120: %o (%o s)\n", g, Realtime() - t;
quit;
