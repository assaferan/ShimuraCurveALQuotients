// Siegel's genus-average identity for the definite binary lattice L_- of X_0^15(2), at a fundamental
// and at conductor discriminants.  With Lpos the positive definite form -<,> on L_-, genus reps L',
//    a(M) = [sum_L' r_L'(M)/|O(L')|] / [sum_L' 1/|O(L')|]      (r counts x with Q^+(x) = M)
// against the product of local densities with the field's L-value,
//    P(M) = (4 pi / sqrt(det S)) * prod_{p in S} alpha_p(M) / ( L(1, chi) prod_{p in S} (1 - chi(p)/p) ),
// chi = chi_{d_0}, S = primes of det S and of M, alpha_p the counting density (p^-k #{Q^+ = M mod p^k}).
// The question is whether the Euler correction (1 - chi(p)/p)^-1 at a conductor prime p | f belongs
// in P (as the identity says) or not (as the data-validated kappa^- normalisation has it).
AttachSpec("ShimuraQuotients.spec");
SetColumns(0);
D := 15; N := 2;
Ld := ShimuraCurveLattice(D, N);
Q := ChangeRing(Ld`Q, Integers()); Qr := ChangeRing(Q, Rationals());
RR := RealField(30);
function density(gram, M, p, kmax)
    // alpha_p(M) for the positive form with <,>-Gram gram: limit of p^-k #{x mod p^k : Q^+(x) = M mod p^k}
    A11 := gram[1,1] div 2; A12 := gram[1,2]; A22 := gram[2,2] div 2;
    last := 0;
    for k in [1..kmax] do
        cnt := 0; md := p^k;
        for u in [0..md-1], t in [0..md-1] do
            if (A11*u^2 + A12*u*t + A22*t^2 - M) mod md eq 0 then cnt +:= 1; end if;
        end for;
        last := cnt / p^k;
    end for;
    return last;
end function;
for d in [-15, -60, -240, -48] do
    lam := ChangeRing(ElementOfNorm(Ld`Q, -d, Ld`O, Ld`basis_L), Integers());
    Lminus := Kernel(Transpose(Matrix(lam*Q)));
    B := BasisMatrix(Lminus); S := B*Q*Transpose(B);           // <,>-Gram, negative definite
    Lpos := LatticeWithGram(-S);
    reps := Representatives(Genus(Lpos));
    mass := &+[RR | 1/#AutomorphismGroup(Lp) : Lp in reps];
    d0 := FundamentalDiscriminant(d); chi := KroneckerCharacter(d0);
    K := QuadraticField(d0); L1 := 2*Pi(RR)*ClassNumber(K)/(#TorsionSubgroup(UnitGroup(K))*Sqrt(RR!-d0));
    detS := Integers()!Determinant(-S);
    printf "\n=== d = %o (d0 %o)  Gram %o  det %o  genus classes %o  mass %o\n", d, d0, Eltseq(-S), Factorization(detS), #reps, mass;
    for M in [1..8] do
        a := &+[RR | Coefficient(ThetaSeries(Lp, 2*M), 2*M)/#AutomorphismGroup(Lp) : Lp in reps] / mass;
        Sset := Sort(SetToSequence(Set(PrimeDivisors(detS)) join Set(PrimeDivisors(M))));
        prod_all := &*[RR | density(-S, M, p, p eq 2 select 9 else 4) : p in Sset];
        corr_all := &*[RR | 1 - Evaluate(chi, p)/p : p in Sset];
        corr_nof := &*[RR | 1 - Evaluate(chi, p)/p : p in Sset | (d div d0) mod p ne 0];
        P_all := 4*Pi(RR)/Sqrt(RR!detS) * prod_all / (L1 * corr_all);
        P_nof := 4*Pi(RR)/Sqrt(RR!detS) * prod_all / (L1 * corr_nof);
        if P_all eq 0 then printf "  M = %o: a(M) = %o   not represented locally\n", M, a; continue; end if;
        printf "  M = %o: a(M) = %o   a/P(all corrections) = %o   a/P(none at p | f) = %o\n", M, a, a/P_all, a/P_nof;
    end for;
end for;
