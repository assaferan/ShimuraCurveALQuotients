SetColumns(0);
PQ<x> := PolynomialRing(Rationals()); PT<T> := PolynomialRing(Rationals());
C := HyperellipticCurve(PQ![9,-33,9,29,-12,6], PQ![0,1,1]);
MS := ModularSymbols(2046, 2); NS := NewSubspace(CuspidalSubspace(MS));
G := [G : G in NewformDecomposition(NS) | Dimension(G) eq 4 and [ Integers() ! (Eltseq(AtkinLehnerOperator(G, p))[1]) : p in [2,3,11,31] ] eq [-1,-1,1,1]][1];
ps := [p : p in PrimesUpTo(30000) | p mod 2046 eq 1];
printf "primes p = 1 mod 2046 below 30000: %o\n", ps;
f := qEigenform(G, ps[#ps] + 1); K := BaseRing(Parent(f));
gcdred := 0;
for p in ps do
    ap := Coefficient(f, p);
    nm := Integers() ! Norm(ap - 1 - p);
    gcdred := GCD(gcdred, nm);
    Lc := PT ! LPolynomial(ChangeRing(C, GF(p)));
    tr := Rationals()!Trace(ap); n := Rationals()!Norm(ap);
    printf "  p = %o: a_p = %o, Norm(a_p - 1 - p) = %o, Euler factor of C matches the orbit: %o\n", p, ap, Factorisation(nm), Lc eq 1 - tr*T + (2*p+n)*T^2 - tr*p*T^3 + p^2*T^4;
end for;
printf "reducibility gcd: %o = %o\n", gcdred, Factorisation(gcdred);
quit;
