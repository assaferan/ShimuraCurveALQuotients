// The Hecke order Z[a_p] of the orbit of the 6_341 star: is it the maximal order Z[sqrt 2]?
// (If so, the star's Jacobian is an O_K-lattice in the Tate module, and the tree argument applies.)
SetColumns(0);
MS := ModularSymbols(2046, 2); NS := NewSubspace(CuspidalSubspace(MS));
G := [G : G in NewformDecomposition(NS) | Dimension(G) eq 4 and [ Integers() ! (Eltseq(AtkinLehnerOperator(G, p))[1]) : p in [2,3,11,31] ] eq [-1,-1,1,1]][1];
f := qEigenform(G, 200); K := BaseRing(Parent(f)); OK := MaximalOrder(K);
r := Roots(PolynomialRing(K)![-2, 0, 1])[1][1];   // sqrt 2 in K
printf "K = %o, sqrt2 = %o\n", DefiningPolynomial(K), r;
coords := [];
for p in [q : q in PrimesUpTo(200) | 2046 mod q ne 0] do
    ap := Coefficient(f, p);
    // write a_p = u + v sqrt2
    v := (ap - Trace(ap)/2) / r; u := Trace(ap)/2;
    assert ap eq u + v*r;
    Append(~coords, <p, u, v>);
end for;
printf "a_p = u + v sqrt2: %o\n", coords[1..12];
printf "gcd of the v's (index of Z[a_p] in Z[sqrt2] over these p): %o\n", GCD([Integers() | c[3] : c in coords]);
printf "a_p with odd v: %o\n", [c : c in coords | IsOdd(Integers()!c[3])][1..Minimum(5, #[c : c in coords | IsOdd(Integers()!c[3])])];
quit;
