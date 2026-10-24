// STEP 5 AT THE RAMIFIED PRIME.  (sqrt 5)^2 = (5), so J[5] is an extension of rho = J[sqrt 5] by
// itself.  If that extension SPLIT as a Galois module, J[5] = rho + rho would contain Galois-stable
// subgroups that are not O_K-submodules -- the kernel of a 25-isogeny to a surface the two gcds
// cannot rule out.  A split J[5] forces dim_{F_5} J(F_p)[5] = 2 dim rho^{Frob_p} to be EVEN at
// every good p; so a single prime with ODD 5-rank proves the extension is not split, after which
// the only Galois-stable subgroups of J[5^n] are the J[sqrt5^k], and every Q-isogeny out of J is
// an isomorphism.  Here J = Jac(C), whose torsion is computable from the curve.
SetColumns(0);
PQ<x> := PolynomialRing(Rationals());
C := HyperellipticCurve(16*x^6+4*x^5+180*x^4-456*x^3+196*x^2+64*x+4, x^2+x);
printf "%-5o %-28o %-10o %o\n", "p", "J(F_p)", "5-rank", "parity";
odd := [];
for p in PrimesInInterval(3, 200) do
    if 2002 mod p eq 0 then continue; end if;
    J := Jacobian(ChangeRing(C, GF(p)));
    G := AbelianGroup(J);
    inv := Invariants(G);
    r5 := #[ i : i in inv | i mod 5 eq 0 ];
    if r5 gt 0 then
        printf "%-5o %-28o %-10o %o\n", p, inv, r5, IsOdd(r5) select "ODD" else "even";
        if IsOdd(r5) then Append(~odd, p); end if;
    end if;
end for;
printf "\nprimes with ODD 5-rank: %o\n", odd;
printf "%o\n", IsEmpty(odd)
    select "=> no odd 5-rank below 200; the splitting question is not settled this way"
    else   "=> J[5] is NOT a split extension: Frobenius at these primes has a Jordan block";
quit;
