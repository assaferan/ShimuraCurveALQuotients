// Which Euler corrections the repo's kappa^- actually applies at a conductor discriminant.
// get_kappa_minus_squared (SchoferFormula.m) divides by prod_{p in S} (1 - chi(p)/p) with
// chi = KroneckerCharacter(d) and S = primes of det(L_-) and of m.  Magma's KroneckerCharacter(d)
// is the PRIMITIVE character of conductor |d_0|, so at a conductor prime p | f in S the correction
// is (1 - chi_{d_0}(p)/p), the field's -- i.e. the literal Euler-product identity.  siegel.m shows
// the same correction in the genus-average identity.  Companion to siegel.m; X_0^15(2).
AttachSpec("ShimuraQuotients.spec");
SetColumns(0);
D := 15; N := 2;
Ld := ShimuraCurveLattice(D, N);
Q := ChangeRing(Ld`Q, Integers());
for d in [-15, -60, -240, -28, -48, -12] do
    lam := ChangeRing(ElementOfNorm(Ld`Q, -d, Ld`O, Ld`basis_L), Integers());
    Lminus := Kernel(Transpose(Matrix(lam*Q)));
    B := BasisMatrix(Lminus); Delta := Determinant(B*Q*Transpose(B));
    d0 := FundamentalDiscriminant(d); f := Isqrt(d div d0);
    chi := KroneckerCharacter(d); chi0 := KroneckerCharacter(d0);
    S := PrimeDivisors(Delta);
    printf "d = %-5o d0 = %-4o f = %o  det L_- = %-14o S = %o  Conductor(KroneckerCharacter(d)) = %o\n", d, d0, f, Factorization(Delta), S, Conductor(chi);
    for p in S do
        printf "    p = %o%o  code's factor 1 - chi(p)/p = %o   field's 1 - chi_d0(p)/p = %o   imprimitive would give %o\n",
            p, f mod p eq 0 select " | f" else "    ", 1 - Evaluate(chi, p)/p, 1 - Evaluate(chi0, p)/p, f mod p eq 0 select 1 else 1 - Evaluate(chi0, p)/p;
    end for;
end for;
