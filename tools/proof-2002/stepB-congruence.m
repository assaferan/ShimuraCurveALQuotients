// STEP B, second half.  Irreducibility of each rho_lambda is not quite enough: at a split l the
// two pieces rho_lambda, rho_lambda' of A[l] could be ISOMORPHIC, and at an inert l the 2-dim
// F_{l^2}-representation could be isomorphic to its Frobenius conjugate; in either case A[l] has
// Galois-stable subgroups that are not O_K-submodules, hence Q-isogenies that are not
// endomorphisms.  Both say a_p = sigma(a_p) mod lambda for all p, i.e. lambda divides the
// congruence modulus between f and its conjugate: lambda | gcd_p (a_p - sigma a_p)^2 = gcd_p (t^2 - 4n).
SetColumns(0);
MS := ModularSymbols(2002, 2); NS := NewSubspace(CuspidalSubspace(MS));
F := 0;
for G in NewformDecomposition(NS) do
    if Dimension(G) ne 4 then continue; end if;
    if [ Integers() ! (Eltseq(AtkinLehnerOperator(G, p))[1]) : p in [2,7,11,13] ] eq [-1,-1,-1,-1]
        then F := G; break; end if;
end for;
g := 0; shown := 0;
for p in PrimesInInterval(3, 400) do
    if 2002 mod p eq 0 then continue; end if;
    ok, q := IsSquare(CharacteristicPolynomial(HeckeOperator(F, p))); assert ok;
    t := -Coefficient(q, 1); n := Coefficient(q, 0);
    d := Integers() ! (t^2 - 4*n);            // (a_p - sigma a_p)^2 = 5 c_p^2
    g := GCD(g, d);
    if shown lt 6 then printf "p = %-4o  t = %-5o n = %-6o (a_p - sigma a_p)^2 = %o\n", p, t, n, d; shown +:= 1; end if;
end for;
printf "...\ngcd over good primes below 400 of (a_p - sigma a_p)^2: %o = %o\n", g, Factorisation(g);
printf "the factor 5 is forced (a_p - sigma a_p is a multiple of sqrt 5); so the congruence\n";
printf "modulus between f and its conjugate is supported on %o\n", [ f[1] : f in Factorisation(g) | f[1] ne 5 ];
// and the one prime that remains, 5 itself: lambda = (sqrt 5) is ramified, rho_lambda is a
// 2-dim F_5-representation and the question is only its irreducibility, settled above by gcd = 1
quit;
