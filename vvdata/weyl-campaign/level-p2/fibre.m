// Which integral-norm cosets of L_-^v/L_- lie in L^v at all, and over which coset of L^v/L?
// Only mu in L^v enter the decomposition of a phi_eta; a coset of the plane outside L^v enters
// no eta, and one lying over eta = 0 is weighted by c_0(0) = 0.  X_0^15(2), the five planes.
AttachSpec("ShimuraQuotients.spec");
D := 15; N := 2;
Ld := ShimuraCurveLattice(D, N);
Q := ChangeRing(Ld`Q, Rationals());           // Gram of the bilinear form on L in basis_L coordinates
n := Nrows(Q);
e := func< i | Vector(Rationals(), [j eq i select 1 else 0 : j in [1..n]]) >;
Lmat := IdentityMatrix(Rationals(), n);        // L = Z^n in these coordinates
for d in [-15, -60, -240, -12, -48] do
    lam := ChangeRing(ElementOfNorm(Ld`Q, -d, Ld`O, Ld`basis_L), Rationals());
    Mx := Matrix(Integers(), n, 1, [Integers() | (e(i)*Q, lam) : i in [1..n]]);
    K := ChangeRing(KernelMatrix(Mx), Rationals());        // rows: basis of L_- in L-coordinates
    gram := K*Q*Transpose(K);
    G := ChangeRing(gram, Integers()); det := Determinant(G); a := Valuation(det, 2); dd := 2^a;
    d0 := FundamentalDiscriminant(d);
    printf "\n=== d = %o (conductor %o, 2 %o): L_N = L_+ + L_- direct? index [L : L_+ + L_-] at 2 ...\n",
        d, Integers()!Sqrt(d/d0), KroneckerSymbol(d0,2) eq 1 select "split" else "INERT";
    // index of L_+ (+) L_- in L: L_+ = Z lam; lattice spanned by lam and K rows
    Lsub := VerticalJoin(Matrix(Rationals(), 1, n, Eltseq(lam)), K);
    idx := AbsoluteValue(Determinant(Lsub)) ;     // vs det of L = 1
    printf "    [L : L_+ (+) L_-] = %o  (2-part 2^%o)\n", idx, Valuation(Integers()!idx, 2);
    seen := {};
    for w in CartesianPower([0..dd-1], 2) do
        v := Vector(Rationals(), [w[1]/dd, w[2]/dd]);
        if not forall{c : c in Eltseq(v*gram) | IsIntegral(c)} then continue; end if;
        key := [c - Floor(c) : c in Eltseq(v)];
        if key in seen or key eq [0,0] then continue; end if;
        Include(~seen, key);
        r := (v*gram, v)/2;
        if r ne Floor(r) then continue; end if;                 // only integral-norm cosets
        mu := Vector(Rationals(), key) * K;                      // the vector in L-coordinates
        inLdual := forall{i : i in [1..n] | IsIntegral((mu*Q, e(i)))};
        inL := forall{c : c in Eltseq(mu) | IsIntegral(c)};
        // which coset of L^v/L: compare Q(mu) mod Z and mu mod L with the nonzero isotropic cosets
        printf "    coset %-12o Q = %-6o  in L^v: %-5o  in L: %-5o  %o\n", key, r, inLdual, inL,
            inLdual select (inL select "-> lies over eta = 0 (weight c_0(0) = 0)" else "-> lies over a NONZERO coset of L^v/L") else "-> enters NO eta";
    end for;
end for;
