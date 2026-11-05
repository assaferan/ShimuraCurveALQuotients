AttachSpec("ShimuraQuotients.spec");
D := 4130; N := 1;
W := { m : m in Divisors(D*N) | GCD(m, (D*N) div m) eq 1 };
printf "genus of the star %o\n", GenusShimuraCurveQuotient(D, N, W);
Y := CreateShimuraQuot(D, N, W); Y`g := GenusShimuraCurveQuotient(D, N, W); Y`CurveID := 0;
curves := GetQuotientsAndGenera([Y]); _ := exists(star){c : c in curves | IsStarCurve(c)};
discs := [-35, -532, -163, -595, -315, -520, -340, -420];
printf "degrees: %o\n", [ <d, DegreeOfFieldOfDefinitionOfCMPoint(star, d)> : d in discs ];
g := 0; t := Realtime();
for p in [q : q in PrimesInInterval(3, 120) | D mod q ne 0] do
    s := [ p^i + 1 - ComputePointsViaTrace(Y, p, i) : i in [1..3] ];
    e := [Rationals() | 1];
    for k in [1..3] do Append(~e, (&+[ (-1)^(i-1) * e[k-i+1] * s[i] : i in [1..k] ]) / k); end for;
    L1 := Integers()!(1 - e[2] + e[3] - e[4] + p*e[3] - p^2*e[2] + p^3);
    g := GCD(g, L1);
    if p le 13 then printf "p = %o: c1 = %o, c2 = %o, c3 = %o, #J(F_p) = %o, T(p^2) char poly x^3 + %o x^2 + %o x + %o\n", p, -e[2], e[3], -e[4], L1, -e[2], e[3] - 3*p, -e[4] + 2*p*e[2]; end if;
end for;
printf "torsion bound gcd #J(F_p), p < 120: %o (%o s)\n", g, Realtime() - t;
quit;
