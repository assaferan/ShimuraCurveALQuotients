AttachSpec("ShimuraQuotients.spec"); SetColumns(0);
P<x> := PolynomialRing(Rationals());
F := 272*x^8 - 72*x^6 - 55*x^4 - 18*x^2 + 17;
H := HyperellipticCurve(F);
X := CreateShimuraQuot(1326, 1, {Integers() | 1,17,78,1326}); X`g := 3;
printf "discriminant: %o\n", Factorisation(Integers()!Discriminant(F));
bad := 0; n := 0;
for p in PrimesInInterval(5, 100) do
    if 1326 mod p eq 0 then continue; end if;
    n +:= 1;
    if #Points(ChangeRing(H, GF(p))) ne ComputePointsViaTrace(X, p, 1) then bad +:= 1; printf "MISMATCH at p = %o\n", p; end if;
end for;
printf "point counts vs trace formula: %o primes, %o mismatches\n", n, bad;
// F_{p^2} counts at a few primes (the full L-polynomial, not only the trace)
for p in [5, 7, 11, 23] do
    n2 := #Points(ChangeRing(H, GF(p^2)));
    printf "p = %o: #H(F_p^2) = %o, trace formula %o\n", p, n2, ComputePointsViaTrace(X, p, 2);
end for;
A, mA := AutomorphismGroup(H); printf "automorphism group order %o\n", #A;
for a in A do
    if Order(a) ne 2 then continue; end if;
    phi := mA(a); eqs := DefiningPolynomials(phi); R := Parent(eqs[1]);
    if eqs[1] eq R.1 and eqs[3] eq R.3 then continue; end if;
    C, q := CurveQuotient(AutomorphismGroup(H, [phi]));
    if Genus(C) ne 1 then continue; end if;
    J := Jacobian(GenusOneModel(C));
    // fixed points of the involution on H: x with phi(x) = x, and the points above them
    num := eqs[1]; den := eqs[3];
    a1 := MonomialCoefficient(num, R.1); b1 := MonomialCoefficient(num, R.3); c1 := MonomialCoefficient(den, R.1); d1 := MonomialCoefficient(den, R.3);
    fx := c1*x^2 + (d1 - a1)*x - b1;   // fixed x-coordinates (x = infinity when c1 = 0)
    flds := [* *];
    for fac in Factorisation(fx) do
        K := NumberField(fac[1]); r := K.1;
        v := Evaluate(F, r); sq := IsSquare(v);
        if sq then Append(~flds, <fac[1], Discriminant(Integers(K))>); else L := ext<K | PolynomialRing(K).1^2 - v>; Append(~flds, <fac[1], Discriminant(Integers(AbsoluteField(L)))>); end if;
    end for;
    if c1 eq 0 then  // infinity is fixed
        lc := LeadingCoefficient(F); Append(~flds, <"oo", FundamentalDiscriminant(SquarefreeFactorization(Integers()!lc))>);
    end if;
    printf "involution x -> (%o x + %o)/(%o x + %o): quotient Jacobian %o; fixed points over fields %o\n", a1, b1, c1, d1, CremonaReference(MinimalModel(J)), flds;
end for;
printf "rational points of small height: %o\n", #Points(H : Bound := 10000);
printf "DONE\n";
