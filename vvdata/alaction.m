// Guo-Yang Thm B assumes O^+_{L,F} = O^+_L: F's components must be constant on the orbits of the
// orthogonal group acting on L^v/L.  There is a short argument that ANY F of the shape
// F_f = sum_gamma (f|gamma) rho(gamma^{-1}) e_0 satisfies this automatically:
//   * sigma in O^+_L is an isometry of the discriminant form, so it COMMUTES with rho_L
//     (rho(T) depends only on Q(eta), rho(S) only on the pairing -- both sigma-invariant);
//   * sigma fixes e_0, since sigma(0) = 0.
// Hence sigma(F_f) = sum_gamma (f|gamma) sigma rho(gamma^{-1}) e_0
//                  = sum_gamma (f|gamma) rho(gamma^{-1}) sigma(e_0) = F_f.
// This script verifies the step the argument turns on -- that the induced action really is an
// isometry of the discriminant form and really does commute with rho_L(S) and rho_L(T) -- using the
// Atkin-Lehner involutions, which generate the relevant group for a Shimura curve
// (O^+_L corresponds to N_{B^x}(O)/Q^x).  Conjugation by a trace-zero w with nrd(w) = m | DN.
AttachSpec("ShimuraQuotients.spec");

for DN in [[6,1],[10,1],[15,2]] do
    D := DN[1]; N := DN[2];
    Ld := ShimuraCurveLattice(D,N);
    Q := ChangeRing(Ld`Q, Rationals()); dn := Ld`denom;
    bas := Ld`basis_L;
    printf "\n=== X0^%o(%o), |G| = %o ===\n", D, N, #Ld`disc_grp;
    CC := ComplexField(40);
    S, Tdiag, elts, i0 := WeilRepresentationComplex(Ld, CC : Dual := true);
    n := #elts;
    idx := AssociativeArray(); for i->g in elts do idx[g] := i; end for;

    for m in [x : x in Divisors(D*N) | x gt 1] do
        ok := true; w := 0;
        try w := ElementOfNorm(Ld`Q, m, Ld`O, bas); catch e ok := false; end try;
        if not ok then printf "  m = %-4o : no trace-zero element of norm %o\n", m, m; continue; end if;
        // conjugation by w on the trace-zero space, in the basis of L
        wq := &+[Eltseq(w)[i]*bas[i] : i in [1..3]];
        rows := [];
        for i in [1..3] do
            y := wq * bas[i] * wq^(-1);
            cs := Solution(Matrix([Eltseq(b) : b in bas]), Vector(Eltseq(y)));
            Append(~rows, Eltseq(cs));
        end for;
        Mat := Matrix(Rationals(), 3, 3, &cat rows);
        isom := (Mat*Q*Transpose(Mat) eq Q);
        integral := &and[IsIntegral(x) : x in Eltseq(Mat)];
        if not (isom and integral) then
            printf "  m = %-4o : conj matrix isometry=%o integral=%o  (skipping)\n", m, isom, integral;
            continue;
        end if;
        // induced permutation of L^v/L
        perm := [];
        for i in [1..n] do
            v := ChangeRing(elts[i]@@Ld`to_disc, Rationals()) * Transpose(Mat);
            Append(~perm, idx[Ld`to_disc(ChangeRing(v, Integers()))]);
        end for;
        isperm := #Set(perm) eq n;
        // does it preserve the quadratic form on the discriminant group?
        vs := [ChangeRing(g@@Ld`to_disc, Rationals()) : g in elts];
        nm := [ (vs[i]*Q, vs[i])/(2*dn^2) : i in [1..n] ];
        keepsQ := &and[ IsIntegral(nm[perm[i]] - nm[i]) : i in [1..n] ];
        // and does it commute with rho?
        P := ZeroMatrix(CC, n, n);
        for i in [1..n] do P[perm[i]][i] := 1; end for;
        commT := Maximum([Abs(Tdiag[perm[i]] - Tdiag[i]) : i in [1..n]]);
        commS := Maximum([Abs(x) : x in Eltseq(P*S - S*P)]);
        printf "  m = %-4o : permutation=%o  preserves Q=%o  [rho(T),sigma]=%o  [rho(S),sigma]=%o  sigma(e_0)=e_%o\n",
               m, isperm, keepsQ, ChangePrecision(commT,4), ChangePrecision(commS,4), perm[i0];
    end for;
end for;
quit;
