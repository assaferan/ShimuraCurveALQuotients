// Validate the general-n D-new trace against modular symbols: the trace of T_n composed with
// W_w on the cuspidal subspace of level DN new at every p | D, averaged over W with the
// Jacquet-Langlands sign.  Grid includes n sharing factors with D, with N, and with both.
AttachSpec("ShimuraQuotients.spec");
function DNewTraceModSym(D, N, n, W)
    M := ModularSymbols(D*N, 2, 1); C := CuspidalSubspace(M);
    for p in PrimeDivisors(D) do C := NewSubspace(C, p); end for;
    if Dimension(C) eq 0 then return 0; end if;
    B := Matrix(Basis(VectorSpace(C))); T := HeckeOperator(M, n);
    tot := 0;
    for w in W do tot +:= (-1)^#PrimeDivisors(GCD(w, D)) * Trace(Solution(B, B*(T * AtkinLehner(M, w)))); end for;
    return tot / #W;
end function;
function ALs(M) return [d : d in Divisors(M) | GCD(d, M div d) eq 1]; end function;
bad := 0; tot := 0; badcop := 0; totcop := 0;
levels := [<1, 12>, <1, 15>, <1, 18>, <1, 20>, <1, 24>, <1, 30>, <1, 36>, <1, 45>, <1, 60>,
           <6, 5>, <6, 7>, <6, 11>, <6, 13>, <10, 3>, <10, 7>, <10, 9>, <14, 3>, <14, 5>, <15, 2>, <15, 4>, <21, 2>, <22, 3>, <26, 3>];
for c in levels do
    D, N := Explode(c); M := D*N;
    als := ALs(M);
    Ws := [{Integers() | 1}, {Integers() | 1, M}] cat [{Integers() | 1, a} : a in als | a ne 1 and a ne M];
    if #als ge 4 then Append(~Ws, Set(als)); end if;
    for W in Ws do
        for n in [2, 3, 4, 5, 6, 8, 9, 10, 12, 25] do
            v := DNewTraceModSym(D, N, n, W);
            f := TraceDNewALFixed(D, N, 2, n, W);
            if GCD(n, M) eq 1 then totcop +:= 1; if v ne f then badcop +:= 1; end if; else tot +:= 1; if v ne f then bad +:= 1; end if; end if;
            if v ne f then printf "BAD D=%o N=%o n=%o W=%o : formula %o modsym %o (gcd %o)\n", D, N, n, W, f, v, GCD(n, M); end if;
        end for;
    end for;
end for;
printf "gcd(n,DN)>1: %o of %o differ;  coprime: %o of %o differ\n", bad, tot, badcop, totcop;
// the cases from the review and the old pins
for c in [<1, 30, 4, {Integers()|1}, -2>, <1, 1848, 7, {Integers()|1,8,231,1848}, -9>, <1, 1848, 49, {Integers()|1,8,231,1848}, -73>] do
    D, N, n, W, expect := Explode(c);
    printf "REVIEW (%o,%o,n=%o): formula %o, modular symbols %o\n", D, N, n, TraceDNewALFixed(D, N, 2, n, W), expect;
end for;
exit;
