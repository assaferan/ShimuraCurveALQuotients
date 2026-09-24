// Structural regression test for the m = 0 term of Schofer's formula.
//
// Schofer ("Borcherds forms and generalizations of singular moduli", Univ. of Maryland thesis 2005 /
// J. reine angew. Math. 629 (2009) 1-36), Theorem 3.3, the (n,2)-Theorem:
//
//     kappa_{eta,lambda}(m_1) = sum_{mu_bar in L/(L_+ + L_-)} sum_{x in eta_+ + mu_+ + L_+}
//                                   kappa^-_{eta_- + mu_-}(m_1 - Q(x)),
//     kappa^-_mu(0) = k_0(0) * phi_{mu}(0),   i.e. ZERO unless mu = 0.
//
// At m_1 = 0 this collapses, for two independent reasons:
//   (1) Q(x) >= 0, so only Q(x) = 0 contributes.  Here L_+ = L cap Q.lambda is RANK ONE and POSITIVE
//       DEFINITE -- for x = c*lambda one has Q(x) = c^2 * (-d) with -d > 0 -- so Q(x) = 0 forces x = 0.
//   (2) phi_mu(0) is a characteristic function, so only eta_- + mu_- = 0 survives.
//
// THE INVARIANT THIS TEST PINS.  Summed over all eta in L^v/L, exactly ONE pair (eta, mu_bar) survives
// both conditions -- and the surviving eta is the TRIVIAL coset.  Consequently:
//   * the m = 0 sum collapses to the single term c_0(0) * k_0(0), and
//   * **the multiplier does not depend on the CM discriminant d**, even though L_+ and L_- very much do.
// The second point is a non-obvious cancellation: the index |L/(L_+ + L_-)| genuinely appears in the
// formula, yet the surviving count is always exactly one.
// The "only the trivial coset" half independently reproduces the delta_{0,mu} that appears in Kudla
// ("Integrals of Borcherds forms", Compositio 137 (2003), Thm II: kappa_mu(0) = (1/2) C delta_{0,mu}),
// in Schofer Thm 3.3 itself, and in Yifan Yang (arXiv:1503.07971, Thm B, where c_0(0) is the WEIGHT of
// the Borcherds form).  Three independent statements of the same fact, now executable.
//
// WHY IT IS WORTH LOCKING IN.  The pipeline's m = 0 multiplier is computed once per form, with no
// dependence on d.  This test is the justification for that design.  It was previously an open worry
// -- a "d-dependent |L/(L_+ + L_-)| factor" was suspected of being a missing ingredient -- and a
// direct computation showed the suspicion was unfounded.  If a future change makes the multiplier
// d-dependent, or breaks the splitting, this test should fail loudly rather than silently shifting
// every CM value.
//
// NOTE ON SCOPE.  This test deliberately checks only the COUNT, which is independent of any
// normalisation convention.  It says nothing about the VALUE of the resulting constant term c_0(0);
// the relative normalisation of the two cusp expansions (Guo-Yang Lemma 24 carries a factor
// M e^{2 pi i/8}/sqrt|L^v/L| on the f|S part) is a separate, still-unsettled question.

// Total m=0 support: #{ (eta, mu_bar) : 0 in eta_+ + mu_+ + L_+  and  eta_- + mu_- = 0 }, together
// with the list of eta that occur.  A plain double loop is used deliberately: |L/(L_+ + L_-)| turns out
// to be tiny (1 on these bases), so this is linear in practice, and an AssociativeArray keyed by
// rational vectors -- the obvious "optimisation" -- is pathologically slow in Magma.
function m0_support(disc_grp, lambda_v, Qint, Qrat, to_disc, denom)
    c_Lplus := Content(lambda_v);
    Lplus := RSpaceWithBasis(Matrix(lambda_v div c_Lplus));
    Lminus := Kernel(Transpose(Matrix(lambda_v*Qint)));
    L := RSpaceWithBasis(IdentityMatrix(Integers(), 3));
    L_quo, L_quo_map := L / (Lplus + Lminus);
    lambda_rat := ChangeRing(lambda_v, Rationals());
    nrm := (lambda_rat*Qrat, lambda_rat);
    total := 0; support := [];
    for eta in disc_grp do
        gamma := ChangeRing(eta@@to_disc, Rationals())/denom;
        c_gp := (gamma*Qrat, lambda_rat)/nrm;
        gamma_minus := gamma - c_gp*lambda_rat;
        for mu_bar in L_quo do
            mu := ChangeRing(mu_bar@@L_quo_map, Rationals());
            c_mp := (mu*Qrat, lambda_rat)/nrm;
            // x = c*lambda with c = c_gp + c_mp + k/c_Lplus; Q(x) = 0 iff c = 0, attainable iff
            // (c_gp + c_mp)*c_Lplus is an integer
            if not IsIntegral((c_gp + c_mp)*c_Lplus) then continue; end if;
            if not IsZero(gamma_minus + (mu - c_mp*lambda_rat)) then continue; end if;
            total +:= 1;
            Append(~support, eta);
        end for;
    end for;
    return total, support;
end function;

procedure test_SchoferM0Structure()
    printf "Testing the m=0 structure of Schofer's formula (Thm 3.3)...";

    // (D, N, [CM discriminants]).  Discriminants are ones the pipeline actually uses on these bases,
    // chosen to span both the "firing" and "non-firing" cases at the level prime.
    // Bases are chosen small on purpose: |L^v/L| = 2(DN)^2, so e.g. X0^10(11) has a discriminant group
    // of ~24200 while these have 1800.  The discriminants are ones that ADMIT an optimal embedding into
    // the order -- ElementOfNorm raises on the others, and each such failure costs seconds, so they are
    // listed explicitly rather than probed.
    cases := [ <15, 2, [-7, -15, -40, -52, -60, -120]>,
               <6,  5, [-4, -19, -24, -40, -51, -120]>,
               <10, 3, [-3, -8, -20, -35, -120]> ];

    nchecked := 0;
    for cs in cases do
        D, N, discs := Explode(cs);
        Ld := ShimuraCurveLattice(D, N);
        Qint := ChangeRing(Ld`Q, Integers());
        Qrat := ChangeRing(Ld`Q, Rationals());
        ndiscs := 0;
        for d in discs do
            lambda := ElementOfNorm(Ld`Q, -d, Ld`O, Ld`basis_L);
            ndiscs +:= 1;

            // (1) L_+ is rank one and positive definite -- the reason Q(x) = 0 forces x = 0.
            lambda_rat := ChangeRing(lambda, Rationals());
            error if (lambda_rat*Qrat, lambda_rat) le 0,
                Sprintf("X0^%o(%o), d=%o: <lambda,lambda> is not positive", D, N, d);

            // (2) summed over all eta, exactly one (eta, mu_bar) pair survives, and its eta is 0
            total, support := m0_support(Ld`disc_grp, lambda, Qint, Qrat, Ld`to_disc, Ld`denom);
            nchecked +:= #Ld`disc_grp;
            error if total ne 1,
                Sprintf("X0^%o(%o), d=%o: total m=0 support count is %o, expected exactly 1 (support %o)",
                        D, N, d, total, support);
            error if not IsZero(support[1]),
                Sprintf("X0^%o(%o), d=%o: the surviving coset is %o, expected the trivial one",
                        D, N, d, support[1]);
        end for;
        error if ndiscs lt 3,
            Sprintf("X0^%o(%o): only %o usable discriminant(s), test would be near-vacuous", D, N, ndiscs);
    end for;

    printf " OK (%o (eta, disc) checks)\n", nchecked;
end procedure;

test_SchoferM0Structure();
