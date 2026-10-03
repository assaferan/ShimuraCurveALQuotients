// The invariants a lemma on the lattice at a conductor prime must reproduce.  For a base (D,N)
// with level prime p || N and discriminants d = f^2 d_0 with p^k || f, p not dividing d_0:
//   the p-part of the index [L : L_+ (+) L_-], the elementary divisors of the plane L_- at p,
//   whether the zero coset's local series has a pole, and how many integral-norm cosets of the
//   plane lie in L^v (only those enter the formula).
AttachSpec("ShimuraQuotients.spec");
P<X> := PolynomialRing(Rationals());
KMAX := 7;
procedure run(D, N, p, ds)
    Ld := ShimuraCurveLattice(D, N);
    Q := ChangeRing(Ld`Q, Rationals()); n := Nrows(Q);
    e := func< i | Vector(Rationals(), [j eq i select 1 else 0 : j in [1..n]]) >;
    printf "\n##### X_0^%o(%o), level prime p = %o #####\n", D, N, p;
    for d in ds do
        ok := true; lam := 0;
        try lam := ChangeRing(ElementOfNorm(Ld`Q, -d, Ld`O, Ld`basis_L), Rationals()); catch err ok := false; end try;
        d0 := FundamentalDiscriminant(d); f := Integers()!Sqrt(d/d0);
        split := KroneckerSymbol(d0, p);
        if not ok then printf "d = %-6o f = %-3o (%o at p): no optimal embedding\n", d, f, split eq 1 select "split" else (split eq -1 select "inert" else "ramified"); continue; end if;
        Mx := Matrix(Integers(), n, 1, [Integers() | (e(i)*Q, lam) : i in [1..n]]);
        K := ChangeRing(KernelMatrix(Mx), Rationals());
        gram := K*Q*Transpose(K); G := ChangeRing(gram, Integers());
        det := Integers()!Determinant(G);
        Lsub := VerticalJoin(Matrix(Rationals(), 1, n, Eltseq(lam)), K);
        idx := Integers()!AbsoluteValue(Determinant(Lsub));
        ed := ElementaryDivisors(G);
        // zero-coset series at p
        alphas := [];
        for k in [1..KMAX] do
            cnt := 0;
            for u, t in [0..p^k-1] do
                x := Vector(Rationals(), [u, t]);
                q := (x*gram, x)/2;
                if IsIntegral(q) and (Integers()!q) mod p^k eq 0 then cnt +:= 1; end if;
            end for;
            Append(~alphas, cnt / p^k);
        end for;
        pole := alphas[KMAX] gt alphas[KMAX-2];
        // integral-norm cosets of the p-part of L_-^v/L_-, and how many lie in L^v
        a := Valuation(det, p); dd := p^a; seen := {}; nint := 0; nLv := 0; nL := 0;
        for w in CartesianPower([0..dd-1], 2) do
            v := Vector(Rationals(), [w[1]/dd, w[2]/dd]);
            if not forall{c : c in Eltseq(v*gram) | IsIntegral(c)} then continue; end if;
            key := [c - Floor(c) : c in Eltseq(v)];
            if key in seen or key eq [0,0] then continue; end if;
            Include(~seen, key);
            r := (v*gram, v)/2;
            if r ne Floor(r) then continue; end if;
            nint +:= 1;
            mu := Vector(Rationals(), key) * K;
            if forall{i : i in [1..n] | IsIntegral((mu*Q, e(i)))} then nLv +:= 1; end if;
            if forall{c : c in Eltseq(mu) | IsIntegral(c)} then nL +:= 1; end if;
        end for;
        printf "d = %-6o d0 = %-5o f = %-3o %-5o | [L : L_+ (+) L_-] p-part p^%o | L_- elem.div. p-vals %o | zero coset %o  %o | integral cosets %o, in L^v %o, in L %o\n",
            d, d0, f, split eq 1 select "split" else "INERT", Valuation(idx, p),
            [Valuation(x, p) : x in ed], alphas[1..5], pole select "POLE" else "no pole", nint, nLv, nL;
    end for;
end procedure;
run(15, 2, 2, [-15, -60, -240, -960, -3, -12, -48, -192, -7, -28, -112]);
run(10, 3, 3, [-8, -72, -648, -387, -603, -43, -67]);
