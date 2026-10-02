// Normalisers of Gamma_0^D(N) for non-squarefree N: the Atkin-Lehner-Newman picture, checked by
// Riemann-Hurwitz.  h = largest divisor of 24 with h^2 | N.  Conjugation by (h 0; 0 1) takes
// Gamma_0^D(h^2 N') to Gamma^D(h) cap Gamma_0^D(N'), normal in Gamma_0^D(N') with quotient
// SL_2(F_h) (S_3 for h = 2, A_4 for h = 3); for 8 | N the same with N' = N/4 and degree 4 over
// Gamma_0^D(2 N'').  The cover is ramified only over the elliptic points of the base.
AttachSpec("ShimuraQuotients.spec");
function genus_and_elliptic(D, N)   // Gamma_0^D(N), any N
    vol := &*[Integers() | p - 1 : p in PrimeDivisors(D)] * N * &*[Rationals() | 1 + 1/p : p in PrimeDivisors(N)];
    e2 := (N mod 4 eq 0) select 0 else &*[Integers() | 1 - KroneckerSymbol(-4, p) : p in PrimeDivisors(D)]
                                        * &*[Integers() | 1 + KroneckerSymbol(-4, p) : p in PrimeDivisors(N)];
    e3 := (N mod 9 eq 0) select 0 else &*[Integers() | 1 - KroneckerSymbol(-3, p) : p in PrimeDivisors(D)]
                                        * &*[Integers() | 1 + KroneckerSymbol(-3, p) : p in PrimeDivisors(N)];
    g := 1 + vol/12 - e2/4 - e3/3;
    return Integers()!g, e2, e3;
end function;

bases := [[6,25],[6,49],[10,9],[14,9],[15,4],[15,8],[21,4],[22,9],[33,4]];
printf "%-6o %-3o %-5o %-9o %-7o %-10o %-12o %o\n", "base", "h", "g", "quotient", "degree", "group", "RH-genus", "full-normaliser quotient genus";
for b in bases do
    D := b[1]; N := b[2];
    g, _, _ := genus_and_elliptic(D, N);
    h := 2^Minimum(Valuation(N, 2) div 2, 3) * 3^Minimum(Valuation(N, 3) div 2, 1);
    if h eq 1 then
        printf "%-6o %-3o %-5o %-9o %-7o %-10o %-12o %o\n", Sprintf("%o_%o", D, N), h, g, "-", "-", "W only", "-", "no extra elements";
        continue;
    end if;
    if h eq 2 and N mod 8 eq 0 then
        Np := N div 4; deg := 4; grp := "(Z/2)^2";
        _, e2, e3 := genus_and_elliptic(D, Np);
        R := 2*e2;     // over an order-2 point: 2 points of index 2; no order-3 points at even level
    elif h eq 2 then
        Np := N div 4; deg := 6; grp := "S_3";
        _, e2, e3 := genus_and_elliptic(D, Np);
        R := 3*e2 + 4*e3;   // order 2: 3 points of index 2; order 3: 2 points of index 3
    else
        Np := N div 9; deg := 12; grp := "A_4";
        _, e2, e3 := genus_and_elliptic(D, Np);
        R := 6*e2 + 8*e3;   // order 2: 6 points of index 2; order 3: 4 points of index 3
    end if;
    gp, _, _ := genus_and_elliptic(D, Np);
    gRH := (deg*(2*gp - 2) + R + 2) div 2;
    // the quotient by the full normaliser: X_0^D(N') / W_{D N'}
    Xs := CreateShimuraQuot(D, Np, Set(Divisors(D*Np)));
    gstar := GenusShimuraCurveQuotient(D, Np, Xs`W);
    printf "%-6o %-3o %-5o %-9o %-7o %-10o %-12o %o\n", Sprintf("%o_%o", D, N), h, g, Sprintf("X_0^%o(%o)", D, Np), deg, grp, gRH, gstar;
end for;
