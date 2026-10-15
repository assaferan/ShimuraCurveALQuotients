// Characteristic polynomials over F_2 of the permutation action of S_{2g+2} on J[2] of a
// hyperelliptic curve of genus g (even subsets of the Weierstrass points modulo the complement):
// for a permutation with cycle lengths c_i the polynomial is prod (x^{c_i}+1) / (x+1)^2.
R<x> := PolynomialRing(GF(2));
for g in [3, 4] do
    n := 2*g + 2; allowed := {};
    for P in Partitions(n) do
        Include(~allowed, (&*[x^c + 1 : c in P]) div (x+1)^2);
    end for;
    printf "genus %o: %o allowed characteristic polynomials out of %o monic of degree %o with f(0)=1\n", g, #allowed, 2^(2*g-1), 2*g;
    printf "  irreducible degree-%o ones allowed: %o\n", 2*g, #[f : f in allowed | IsIrreducible(f)];
end for;
quit;
