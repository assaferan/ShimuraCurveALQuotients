// Regression test for EisensteinLocalFactors.m -- Kudla-Yang's local Eisenstein factors written out
// from the paper (Sci. China Math. 53 (2010) 2275-2316) rather than read off the lattice.
//
// The decisive check is (2): wherever the local quaternion algebra is SPLIT and the lattice is
// unimodular (p not dividing D*N), Proposition 5.3 must reproduce the pipeline's own
// LocalWhittakerAtOne.  It does, on every base tested.  That pins (2.16), the L-factor, the local
// zeta factor, the conductor indexing and the character convention all at once -- which is what makes
// the RAMIFIED branch (2.17) trustworthy at p | D, where no independent check is available.

procedure test_bp_special_values()
    printf "  b_p closed-form spot checks...";
    // Ramified, conductor exponent k = 0, p ramified in Q(sqrt(kappa m)) so v_p = 0.  Then (2.17)
    // collapses: numerator = (1 - p^2 X^2) + p^2 X^2 - p^2 X^4 = 1 - p^2 X^4, and dividing by
    // 1 - p X^2 leaves 1 + p X^2.
    for p in [3, 5, 7] do
        // kappa m = -p gives d = -p or -4p, ramified at p, and c = 1 so k = 0
        km := Rationals()!(-p);
        if KYSplittingType(km, p) ne 0 then continue; end if;
        for X in [Rationals() | 1/2, 1/3, 1/p] do
            assert KYbp(km, p, true, X) eq 1 + p*X^2;
        end for;
    end for;
    // At a prime not dividing D and not dividing the conductor, (2.16) is identically 1: with k = 0
    // the numerator is 1 - vX + vX - pX^2 = 1 - pX^2, which is exactly the denominator.
    for p in [7, 11, 13] do
        for m in [1..12] do
            km := Rationals()!m;
            if p gt 3 and (4*m) mod p ne 0 then
                assert KYbp(km, p, false, Rationals()!1/p) eq 1;
            end if;
        end for;
    end for;
    printf " ok\n";
end procedure;

procedure test_prop53_matches_pipeline_at_good_primes()
    printf "  Prop 5.3 vs LocalWhittakerAtOne at good primes...";
    Lfull := RSpaceWithBasis(IdentityMatrix(Integers(), 3));
    zero := Vector(Rationals(), [0,0,0]);
    checked := 0;
    for DN in [ <15,2>, <6,5>, <10,3>, <6,7>, <10,11> ] do
        D := DN[1]; N := DN[2];
        Ld := ShimuraCurveLattice(D, N);
        negQ := -ChangeRing(Ld`Q, Integers());
        for m in [1..30] do
            for p in PrimeDivisors(2*3*5*7*11*13*m) do
                if p gt 30 then continue; end if;
                if (D*N) mod p eq 0 then continue; end if;          // good primes only
                assert LocalWhittakerAtOne(Rationals()!m, p, zero, Lfull, negQ)
                       eq KYWhittaker53(Rationals()!m, p, false);
                checked +:= 1;
            end for;
        end for;
    end for;
    assert checked ge 200;
    printf " ok (%o checks)\n", checked;
end procedure;

procedure test_prop54_matches_pipeline_at_the_level_prime()
    printf "  Prop 5.4 vs LocalWhittakerAtOne at p | N...";
    Lfull := RSpaceWithBasis(IdentityMatrix(Integers(), 3));
    zero := Vector(Rationals(), [0,0,0]);
    checked := 0;
    for DN in [ <6,5>, <6,7>, <10,3>, <10,11>, <22,7> ] do
        D := DN[1]; N := DN[2];
        Ld := ShimuraCurveLattice(D, N);
        negQ := -ChangeRing(Ld`Q, Integers());
        for m in [1..30] do
            assert LocalWhittakerAtOne(Rationals()!m, N, zero, Lfull, negQ)
                   eq KYWhittaker54(Rationals()!m, N, 1);
            checked +:= 1;
        end for;
    end for;
    printf " ok (%o checks)\n", checked;
end procedure;

procedure test_ramified_branch_differs()
    printf "  the ramified branch is a genuinely different factor...";
    // The point of (2.17): at p | D the anisotropic density is NOT the split one.  Exhibit the
    // divergence explicitly so that a refactor collapsing the two cases fails here.
    ndiff := 0;
    for p in [2, 3, 5] do
        for m in [1..20] do
            km := Rationals()!m;
            if KYbp(km, p, true, Rationals()!1/p) ne KYbp(km, p, false, Rationals()!1/p) then
                ndiff +:= 1;
            end if;
        end for;
    end for;
    assert ndiff ge 20;
    printf " ok (%o of 60 indices differ)\n", ndiff;
end procedure;

// The local Whittaker of the N-SCALED binary lattice L_- at the level prime, in closed form.
// This is the ingredient a level-N analogue of Schofer's Thm 4.1 needs, and the reason that theorem
// cannot simply be transplanted: at ord_N(m) = 0 the factor is 1 - X = 1 - N^(-s), the INVERSE LOCAL
// ZETA FACTOR, whose zero at s = 0 is a normalisation artifact and not a failure of the space to
// represent m.  Schofer Thm 4.1 / Yang Lemma 16 / kappaminus all treat a vanishing local factor as a
// representability obstruction, which is why they discard these terms.
procedure test_binary_level_whittaker()
    printf "  W_N of the N-scaled binary L_- in closed form...";
    checked := 0;
    for DN in [ <10,3,-8>, <15,2,-7>, <6,5,-19> ] do
        D := DN[1]; N := DN[2]; d := DN[3];
        Ld := ShimuraCurveLattice(D, N);
        Q := ChangeRing(Ld`Q, Integers());
        lam := ElementOfNorm(Ld`Q, -d, Ld`O, Ld`basis_L);
        Lminus := Kernel(Transpose(Matrix(lam*Q)));
        Bm := BasisMatrix(Lminus);
        // at a firing discriminant L_- is N-scaled, not unimodular, at the level prime
        assert Valuation(Determinant(Bm*Q*Transpose(Bm)), N) eq 2;
        zero := Vector(Rationals(), [0,0,0]);
        R<X> := PolynomialRing(Rationals());
        // The closed form, for EVERY valuation j = ord_N(m):
        //     W_N(X) = 1 + (N-1)(X + X^2 + ... + X^j) - X^(j+1)
        // so  W_N(1) = (N-1)*ord_N(m)  -- vanishing exactly when N does not divide m, which IS the
        // empirical level support rule -- and W_N'(1) = (j+1)((N-1)j/2 - 1).
        // At j = 0 this is 1 - X = 1 - N^(-s), the inverse local zeta factor, whose zero at s = 0 is
        // a normalisation artifact rather than a failure to represent m.
        for m in [1..60] do
            j := Valuation(m, N);
            W := R ! LocalWhittakerPolynomial(Rationals()!m, N, zero, Lminus, Q);
            pred := 1 - X^(j+1) + (N-1)*&+[R | X^i : i in [1..j]];
            assert W eq pred;
            assert Evaluate(W, 1) eq (N-1)*j;
            assert Evaluate(W, 1) eq LocalWhittakerAtOne(Rationals()!m, N, zero, Lminus, Q);
            assert Evaluate(Derivative(W), 1) eq (j+1)*((N-1)*j/2 - 1);
            checked +:= 1;
        end for;
    end for;
    assert checked ge 150;
    printf " ok (%o checks)\n", checked;
end procedure;

procedure test_EisensteinLocalFactors()
    printf "Testing Kudla-Yang local Eisenstein factors...\n";
    test_bp_special_values();
    test_prop53_matches_pipeline_at_good_primes();
    test_prop54_matches_pipeline_at_the_level_prime();
    test_ramified_branch_differs();
    test_binary_level_whittaker();
    printf "Done!\n";
end procedure;

test_EisensteinLocalFactors();
