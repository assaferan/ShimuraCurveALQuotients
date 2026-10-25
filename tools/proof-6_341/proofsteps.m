// THE CHEAP STEPS OF THE PR #139 TEMPLATE FOR THE STAR OF X_0^6(341): the newform orbit, the Euler
// factors, the conductor (Magma's, to be proved at 2 by CHIMP), the reducibility and congruence
// gcds, and the parity test at the ramified prime of the RM field.
SetColumns(0);
PQ<x> := PolynomialRing(Rationals()); PT<T> := PolynomialRing(Rationals());
C := HyperellipticCurve(PQ![9,-33,9,29,-12,6], PQ![0,1,1]);
printf "C : %o\n", C;
printf "discriminant: %o\n", Factorisation(Integers() ! Discriminant(C));
cond := Conductor(C);
printf "conductor of Jac(C) (Magma; Ogg at 2 without guarantee): %o = %o   (2046^2 = %o)\n", cond, Factorisation(cond), 2046^2;
good := [p : p in PrimesUpTo(60) | 2046 mod p ne 0];
Lc := [ <p, PT ! LPolynomial(ChangeRing(C, GF(p)))> : p in good ];
// the newform orbit: level 2046, dimension 2, coefficient field Q(sqrt 2), Atkin-Lehner signs
// -1 at the primes of D = 6 (the Jacquet-Langlands sign flip) and +1 at the primes of N = 341
MS := ModularSymbols(2046, 2); NS := NewSubspace(CuspidalSubspace(MS));
printf "dim S_2^new(2046) = %o\n", Dimension(NS);
cands := [];
for G in NewformDecomposition(NS) do
    if Dimension(G) ne 4 then continue; end if;   // 2-dim orbit (dimension 4 as a +/- space)
    signs := [ Integers() ! (Eltseq(AtkinLehnerOperator(G, p))[1]) : p in [2,3,11,31] ];
    f := qEigenform(G, 70); K := BaseRing(Parent(f));
    printf "orbit: AL signs %o at [2,3,11,31], coefficient field %o\n", signs, DefiningPolynomial(K);
    // Euler factors at good primes against Jac(C)
    ok := true;
    for t in Lc do
        p := t[1]; ap := Coefficient(f, p);
        tr := Rationals() ! Trace(ap); nm := Rationals() ! Norm(ap);
        Lf := 1 - tr*T + (2*p + nm)*T^2 - tr*p*T^3 + p^2*T^4;     // the product over the two conjugates
        if Lf ne t[2] then ok := false; end if;
    end for;
    printf "   Euler factors match Jac(C) at %o: %o\n", good, ok;
    if ok then Append(~cands, G); end if;
end for;
printf "%o matching orbit(s)\n", #cands;
if #cands eq 1 then
    G := cands[1]; f := qEigenform(G, 400); K := BaseRing(Parent(f)); OK := MaximalOrder(K);
    printf "K = %o, discriminant %o\n", K, Discriminant(OK);
    // reducibility: lambda | Norm(a_p - 1 - p) for p = 1 mod 2046 (characters of conductor dividing 2046)
    gcdred := 0;
    for p in [q : q in PrimesUpTo(400) | q mod 2046 eq 1] do gcdred := GCD(gcdred, Integers() ! Norm(Coefficient(f, p) - 1 - p)); end for;
    ps := [q : q in PrimesUpTo(400) | q mod 2046 eq 1];
    printf "reducibility gcd over p = 1 mod 2046 (%o): %o\n", ps, gcdred;
    // congruence modulus between f and its conjugate: gcd of (a_p - sigma a_p)^2 = t^2 - 4n
    gcdcong := 0;
    for p in [q : q in PrimesUpTo(400) | 2046 mod q ne 0] do
        ap := Coefficient(f, p); gcdcong := GCD(gcdcong, Integers() ! (Trace(ap)^2 - 4*Norm(ap)));
    end for;
    printf "congruence gcd (t^2 - 4n over good p < 400): %o = %o\n", gcdcong, Factorisation(gcdcong);
    // the ramified prime of K = Q(sqrt 2) is 2: parity of the 2-rank of J(F_p)
    J := Jacobian(C);
    printf "2-ranks of J(F_p): %o\n", [ <p, Valuation(#AbelianGroup(BaseChange(J, GF(p))), 2), #Invariants(AbelianGroup(BaseChange(J, GF(p)))) > : p in [5,7,13,17,19,23] ];
    printf "J(F_p) invariants: %o\n", [ <p, Invariants(AbelianGroup(BaseChange(J, GF(p))))> : p in [5,7,13,17,19,23,29,37] ];
end if;
quit;
