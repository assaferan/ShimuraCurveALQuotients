// Ramification of the degeneracy map X_0^D(N) -> X_0^D(N'), N = h^2 N'.  The map of orbifold
// quotients H/Gamma_0(N) -> H/Gamma_0(N') ramifies exactly where the stabiliser shrinks, i.e. over
// the base's ELLIPTIC points, and an elliptic point of order e on the base has n_top/n_bot elliptic
// preimages (index 1) and (deg - n_top/n_bot)/e further ones of index e.  So
//     R = sum_{e in {2,3}} n_bot(e) * (deg - n_top(e)/n_bot(e)) * (e-1)/e,
// and Riemann-Hurwitz must then reproduce the genus formula's g_top.  No Atkin-Lehner quotient is
// involved, so this tests the law on its own.
AttachSpec("ShimuraQuotients.spec");
psi := func<M | M * &*[Rationals() | 1 + 1/q : q in PrimeDivisors(M)]>;
function ell(D, N)
    e2 := (N mod 4 eq 0) select 0 else &*[Integers() | 1 - KroneckerSymbol(-4, p) : p in PrimeDivisors(D)]
                                        * &*[Integers() | 1 + KroneckerSymbol(-4, p) : p in PrimeDivisors(N)];
    e3 := (N mod 9 eq 0) select 0 else &*[Integers() | 1 - KroneckerSymbol(-3, p) : p in PrimeDivisors(D)]
                                        * &*[Integers() | 1 + KroneckerSymbol(-3, p) : p in PrimeDivisors(N)];
    return e2, e3;
end function;
function gen(D, N)
    vol := &*[Integers() | p - 1 : p in PrimeDivisors(D)] * N * &*[Rationals() | 1 + 1/p : p in PrimeDivisors(N)];
    e2, e3 := ell(D, N);
    return Integers() ! (1 + vol/12 - e2/4 - e3/3);
end function;
printf "%-9o %-6o %-6o %-5o %-18o %-18o %-7o %-7o %o\n",
       "base", "g_top", "g_bot", "deg", "elliptic bot (2,3)", "elliptic top (2,3)", "R", "g pred", "verdict";
for b in [[15,4,1],[21,4,1],[33,4,1],[15,8,2],[10,9,1],[14,9,1],[22,9,1],[6,25,1],[6,49,1],[10,3,1],[14,3,1]] do
    D := b[1]; N := b[2]; Np := b[3];
    deg := Integers() ! (psi(N) / psi(Np));
    gb := gen(D, Np); e2b, e3b := ell(D, Np); e2t, e3t := ell(D, N);
    R := 0;
    for pair in [<2, e2b, e2t>, <3, e3b, e3t>] do
        e := pair[1]; nb := pair[2]; nt := pair[3];
        if nb eq 0 then continue; end if;
        R +:= nb * (deg - nt/nb) * (e - 1)/e;
    end for;
    gpred := (deg*(2*gb - 2) + R + 2) / 2;
    g_top := gen(D, N);
    printf "%-9o %-6o %-6o %-5o %-18o %-18o %-7o %-7o %o\n", Sprintf("%o_%o", D, N), g_top, gb, deg,
           Sprintf("%o, %o", e2b, e3b), Sprintf("%o, %o", e2t, e3t), R, gpred,
           gpred eq g_top select "OK" else "MISMATCH";
end for;
