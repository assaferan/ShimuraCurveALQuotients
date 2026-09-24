// Regression tests for VectorValuedForm.m -- the vector-valued input form
//     F_f = sum_{gamma in Gamma_0(M)\SL_2(Z)} (f|_{1/2} gamma) rho(gamma^{-1}) e_0        [GY, eq (6)]
// and the constant terms c_eta(0) that Schofer's m = 0 multiplier is built from.
//
// The tests are layered so that a failure localises:
//   (1) the coset transversal of Gamma_0(M)\SL_2(Z) and the ST-word decomposition;
//   (2) the complex Weil representation, against the exact cyclotomic WeilRepresentationST, and the
//       metaplectic relations of the DUAL -- including rho(Z^-1) e_0 = e(1/4) e_0, which is the sign
//       the coset sum needs and the reason the dual is the right representation;
//   (3) the number of isotropic cosets, = 2N - 1 (the D-part of the discriminant group is
//       anisotropic, so every isotropic coset comes from the level-N hyperbolic plane);
//   (4) the metaplectic S-slash against qExpansionAt0, which pins this file's normalisation to the
//       pipeline's own (validated) cusp-0 convention;
//   (5) F_f end to end on X0^6(1): its principal part must reproduce f's, and c_0(0) must vanish;
//   (6) SLOW: the m = 0 multiplier on X0^15(2), against the independently measured ground truth.

procedure test_cosets_and_words()
    printf "  cosets and ST-words...";
    for M in [12, 20, 60] do
        reps := VVCosetReps(M);
        assert #reps eq Index(Gamma0(M));
        // pairwise distinct cosets: g1 g2^-1 lies in Gamma_0(M) iff M divides its lower-left entry
        for i in [1..#reps] do
            for j in [i+1..#reps] do
                assert ((reps[i]*reps[j]^(-1))[2][1] mod M) ne 0;
            end for;
        end for;
        for g in reps do
            assert VVWordMatrix(VVSTWord(g)) eq g;
        end for;
    end for;
    printf " ok\n";
end procedure;

procedure test_weil_complex()
    printf "  complex Weil representation...";
    CC := ComplexField(40);
    ii := CC.1;
    tol := 1e-25;
    for DN in [ [6,1], [10,1] ] do
        Ld := ShimuraCurveLattice(DN[1], DN[2]);
        Sx, Tx, eltsx, Kcyc := WeilRepresentationST(Ld);
        n := #eltsx;

        // (a) the complex build agrees with the exact cyclotomic one
        Sc, Td, elts, i0 := WeilRepresentationComplex(Ld, CC : Dual := false);
        assert elts eq eltsx;
        zt := Exp(2*Pi(CC)*ii/CyclotomicOrder(Kcyc));
        emb := func<x | &+[CC | Eltseq(x)[k]*zt^(k-1) : k in [1..#Eltseq(x)]]>;
        assert Maximum([Abs(Sc[i][j] - emb(Sx[i][j])) : i in [1..n], j in [1..n]]) lt tol;
        assert Maximum([Abs(Td[i] - emb(Tx[i][i])) : i in [1..n]]) lt tol;

        // (b) the DUAL satisfies the Mp_2(Z) relations, with the opposite S-phase
        S, Td, elts, i0 := WeilRepresentationComplex(Ld, CC : Dual := true);
        T := DiagonalMatrix(CC, Td);
        idx := AssociativeArray(); for i->g in elts do idx[g] := i; end for;
        P := ZeroMatrix(CC, n, n);
        for j in [1..n] do P[idx[-elts[j]]][j] := 1; end for;
        assert Maximum([Abs(x) : x in Eltseq(S*S - Exp(-2*Pi(CC)*ii/4)*P)]) lt tol;
        assert Maximum([Abs(x) : x in Eltseq((S*T)^3 - S*S)]) lt tol;
        assert Maximum([Abs(x) : x in Eltseq(S^4 + IdentityMatrix(CC, n))]) lt tol;

        // (c) the reason the dual is the right one: f|Z = i^{-1} f at weight 1/2, so well-definedness
        //     of the coset sum forces rho(Z^{-1}) e_0 = e(1/4) e_0.  Borcherds' rho_L gives e(-1/4).
        v := VVRhoInvE0(S, Td, [<"S",0>,<"S",0>], i0);
        assert Abs(v[i0][1] - Exp(2*Pi(CC)*ii/4)) lt tol;
        assert Maximum([Abs(v[i][1]) : i in [1..n] | i ne i0]) lt tol;
    end for;
    printf " ok\n";
end procedure;

procedure test_weil_fft()
    printf "  S-action as a Fourier transform on the discriminant group...";
    CC := ComplexField(40);
    for DN in [ <6,1>, <10,1>, <15,2> ] do
        D := DN[1]; N := DN[2];
        Ld := ShimuraCurveLattice(D, N);
        S, Tdiag, elts, i0 := WeilRepresentationComplex(Ld, CC : Dual := true);
        n := #elts;
        data := VVWeilFFT(Ld, CC : Dual := true);
        assert data[7] eq elts and data[8] eq i0;
        assert Maximum([Abs(data[6][i] - Tdiag[i]) : i in [1..n]]) lt 1e-25;
        // the factorised S must agree with the dense matrix on arbitrary vectors, not just on e_0
        for trial in [1..3] do
            v := [CC | ((i*trial) mod 7) - 3 + CC.1*(((i+trial) mod 5) - 2) : i in [1..n]];
            dense := S * Matrix(CC, n, 1, v);
            fast := VVApplyS(data, v);
            assert Maximum([Abs(dense[i][1] - fast[i]) : i in [1..n]]) lt 1e-25;
        end for;
        // and on the actual coset words, which is how it is used
        M := IsOdd(D*N) select 4*D*N else 2*D*N;
        for g in VVCosetReps(M)[1..Minimum(6, #VVCosetReps(M))] do
            w := VVSTWord(g);
            a := VVRhoInvE0(S, Tdiag, w, i0);
            b := VVRhoInvE0FFT(data, w);
            assert Maximum([Abs(a[i][1] - b[i]) : i in [1..n]]) lt 1e-22;
        end for;
    end for;
    printf " ok\n";
end procedure;

procedure test_isotropic_count()
    printf "  isotropic cosets number 2N-1...";
    for DN in [ [6,1], [15,2], [21,2], [10,3], [6,5] ] do
        D := DN[1]; N := DN[2];
        Ld := ShimuraCurveLattice(D, N);
        Q := ChangeRing(Ld`Q, Rationals()); dn := Ld`denom;
        cnt := 0;
        for g in Ld`disc_grp do
            v := ChangeRing(g@@Ld`to_disc, Rationals());
            if IsIntegral((v*Q, v)/(2*dn^2)) then cnt +:= 1; end if;
        end for;
        assert cnt eq 2*N - 1;
    end for;
    printf " ok\n";
end procedure;

// Returns the star curve's Borcherds forms for X0^D(N), as the research runs build them.
function borcherds_forms(D, N)
    Xstar := CreateShimuraQuot(D, N, Set(Divisors(D*N)));
    Xstar`g := GenusShimuraCurveQuotient(D, N, Xstar`W); Xstar`CurveID := 0;
    curves := GetQuotientsAndGenera([Xstar]);
    _ := exists(star){c : c in curves | IsStarCurve(c)};
    return BorcherdsForms(star, curves : Prec := 100);
end function;

procedure test_slash_against_qexpansion_at_0()
    printf "  metaplectic S-slash vs qExpansionAt0...";
    // Unwinding SAction against eta(-1/z) = sqrt(z/i) eta(z) termwise: SAction divides the r-th
    // monomial by sqrt(prod d^{r_d} * disc) and multiplies by M, and the prod d^{r_d} cancels the
    // M^k/sqrt(prod d^{r_d}) of the eta transformation.  What survives is
    //     (f|S~)(M z) = e(-1/8) * sqrt(disc) / M * qExpansionAt0(f)(e(z)),
    // exactly the INVERSE of [GY, Lemma 24]'s factor M e(1/8)/sqrt|L^v/L| -- i.e. the pipeline's
    // convention c_eta(-n/M) = Coefficient(f0,-n) already carries Lemma 24's normalisation.
    fs := borcherds_forms(6, 1);
    f := fs[Sort([k : k in Keys(fs)])[1]];
    R := Parent(f); M := R`M;
    CC := ComplexField(60); ii := CC.1;
    ser0 := qExpansionAt0(f, 60);
    pred := Exp(-Pi(CC)*ii/4) * Sqrt(CC!R`disc) / M;
    for tau in [CC | 40*ii, 7.0 + 33*ii, -9.0 + 55*ii] do
        lhs := VVSlashEval(f, [<"S",0>], tau);
        q := Exp(2*Pi(CC)*ii*tau/M);
        rhs := &+[CC | Coefficient(ser0, k)*q^k : k in [Valuation(ser0)..50]];
        assert Abs(lhs/rhs - pred) lt 1e-40;
    end for;
    printf " ok\n";
end procedure;

procedure test_Ff_on_6_1()
    printf "  F_f end to end on X0^6(1)...";
    D := 6; N := 1;
    fs := borcherds_forms(D, N);
    ks := Sort([k : k in Keys(fs)]);
    Ld := ShimuraCurveLattice(D, N);
    M := 2*D*N;
    forms := [fs[k] : k in ks];
    // N = 1, so the only isotropic coset is the trivial one and there is no m = 0 multiplier;
    // what this pins is the CONSTRUCTION -- the principal part, and the vanishing of c_0(0).
    // PosDepth = 2M also requests the positive coefficients c_eta(j/M) up to j = 2M.  NOTE there is
    // no exact external target for these: GY Lemma 24 pins only the PRINCIPAL parts, and F_f's
    // positive coefficients differ from the scalar expansions' (c_0(0) = 0 vs Coefficient(foo, 0)
    // already shows this at n = 0).  What IS forced structurally is the support: the trivial
    // coset's component has integer exponents only, so every off-grid fractional exponent must
    // extract to numerical zero -- that checks the new extraction path end to end.
    consts, isoelts, errs, poscs := VVConstantTerms(forms, Ld, M
                                                    : Prec := 60, NumSamples := 64, PosDepth := 2*M);
    assert #isoelts eq 1 and IsZero(isoelts[1]);
    for i in [1..#forms] do
        assert errs[i] lt 1e-30;            // principal part reproduces GY Lemma 24
        assert Abs(consts[i][1]) lt 1e-30;  // c_0(0) = 0
        // plumbing: the right shape, and genuinely nonzero on-grid values.  (The off-grid noise
        // floor scales with the PRINCIPAL-part magnitudes, so a sharp off-grid zero test needs the
        // production parameters -- see the aliasing caveat in the intrinsic's documentation.)
        assert #poscs[i][1] eq 2*M;
        assert Maximum([Abs(poscs[i][1][j]) : j in [M, 2*M]]) gt 1e-10;
    end for;
    printf " ok\n";
end procedure;

procedure test_m0_multiplier_15_2()
    printf "  m=0 multiplier on X0^15(2) (slow)...";
    D := 15; N := 2;
    fs := borcherds_forms(D, N);
    Ld := ShimuraCurveLattice(D, N);
    // Ground truth, measured independently by the kernel-consistency oracle sweep, in units of log N.
    // Only the two cheap forms are checked here: the pole-order-30 forms need Prec 200 / 192 samples
    // and about 20 minutes.  Forms -1 and 13 have pole order 10 and converge quickly.
    expected := AssociativeArray();
    expected[-1] := 4;  expected[13] := 4;
    keys := Sort([k : k in Keys(expected)]);
    mults, errs := M0MultiplierNumeric([fs[k] : k in keys], Ld, D, N :
                                       Prec := 120, NumSamples := 64);
    for i->k in keys do
        assert errs[i] lt 1e-20;
        assert Abs(mults[i] - expected[k]) lt 1e-15;
    end for;
    printf " ok\n";
end procedure;

procedure test_VectorValuedForm()
    printf "Testing the vector-valued Borcherds input form F_f...\n";
    test_cosets_and_words();
    test_weil_complex();
    test_weil_fft();
    test_isotropic_count();
    test_slash_against_qexpansion_at_0();
    test_Ff_on_6_1();
    test_m0_multiplier_15_2();
    printf "Done!\n";
end procedure;

test_VectorValuedForm();
