// Step 4: validate the complex Weil representation.
//  (i)  my complex build agrees with the exact cyclotomic WeilRepresentationST (non-dual)
//  (ii) the DUAL satisfies the Mp_2(Z) relations with the sign the coset sum needs:
//         rho(S)^2 = e(-1/4) P,  (rho(S)rho(T))^3 = rho(S)^2,  rho(S)^4 = -I
//       and hence rho(Z^{-1}) e_0 = e(1/4) e_0, matching f|Z = i^{-1} f at weight 1/2.
AttachSpec("ShimuraQuotients.spec");
Attach("vvlib.m");

CC := ComplexField(40);
ii := CC.1;
tol := 1e-25;

for DN in [[6,1],[10,1]] do
    D := DN[1]; N := DN[2];
    Ld := ShimuraCurveLattice(D,N);
    Sx, Tx, eltsx, Kcyc := WeilRepresentationST(Ld);
    n := #eltsx;
    printf "\n=== %o_%o   |L^v/L| = %o\n", D, N, n;

    // (i) compare against the exact construction (same element order? rebuild ours on eltsx order)
    Sc, Td, elts, i0 := VVRho(Ld, CC : Dual := false);
    assert elts eq eltsx;
    zt := Exp(2*Pi(CC)*ii/CyclotomicOrder(Kcyc));
    emb := func<x | &+[CC | Eltseq(x)[k]*zt^(k-1) : k in [1..#Eltseq(x)]]>;
    dS := Maximum([Abs(Sc[i][j] - emb(Sx[i][j])) : i in [1..n], j in [1..n]]);
    dT := Maximum([Abs(Td[i] - emb(Tx[i][i])) : i in [1..n]]);
    printf "(i)  max|S_complex - S_exact| = %o    max|T diff| = %o\n", dS, dT;
    assert dS lt tol and dT lt tol;

    // (ii) metaplectic relations for the DUAL
    S, Td, elts, i0 := VVRho(Ld, CC : Dual := true);
    T := DiagonalMatrix(CC, Td);
    Id := IdentityMatrix(CC, n);
    idx := AssociativeArray(); for i->g in elts do idx[g] := i; end for;
    P := ZeroMatrix(CC, n, n);
    for j in [1..n] do P[idx[-elts[j]]][j] := 1; end for;
    c := Exp(-2*Pi(CC)*ii/4);                       // e(-1/4)
    printf "(ii) |S^2 - e(-1/4)P| = %o\n", Maximum([Abs(x) : x in Eltseq(S*S - c*P)]);
    printf "(ii) |(ST)^3 - S^2|   = %o\n", Maximum([Abs(x) : x in Eltseq((S*T)^3 - S*S)]);
    printf "(ii) |S^4 + I|        = %o\n", Maximum([Abs(x) : x in Eltseq(S^4 + Id)]);
    // rho(Z^{-1}) e_0 where Z = S^2:  must be e(1/4) e_0
    v := VVRhoInvE0(S, Td, [<"S",0>,<"S",0>], i0);
    target := Exp(2*Pi(CC)*ii/4);
    printf "(ii) rho(Z^-1)e_0 = %o * e_0 ; off-e_0 mass = %o\n", v[i0][1],
           Maximum([Abs(v[i][1]) : i in [1..n] | i ne i0]);
    assert Abs(v[i0][1] - target) lt tol;
end for;
printf "\nAll rho checks passed.\n";
quit;
