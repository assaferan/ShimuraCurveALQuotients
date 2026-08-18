// Regression test for the principled outer-m=0 multiplier sum_eta c_eta(0) of Schofer's formula.
//
// The m=0 level-prime term of Schofer's formula (Yang, arXiv:1503.07971, eq (11)-(12) + Lemma 20) is
// (sum_eta c_eta(0)) * kappa_0^-(0), where sum_eta c_eta(0) is the constant term of the vector-valued
// Borcherds input F_f. M0Multiplier computes it via the Borcherds/Eisenstein obstruction: pair F_f
// (weight 1/2, rho_L) against the holomorphic weight-3/2 Eisenstein series for the dual rho_L^*, whose
// Fourier coefficients are the Kudla-Yang local Whittaker products (Sci China Math 53 (2010), Prop
// 2.6(ii) + Prop 5.3) times the class-number/L-value normalization, all collapsing to -96/(D*N)*T.
//
// This pins the coefficient's exact normalization -- LocalWhittakerAtOne, the chi_{-m} character, the
// class number h/w, and the -96/(D*N) prefactor -- independently of the full Borcherds/CM pipeline
// (which SchoferIsometry.m exercises end-to-end and is much slower).
//
// KNOWN LIMITATION (see the commit that introduced m0_multiplier): the coefficient is exact only for
// single-surviving-term inputs. It is validated here on X0^15(2), whose star hauptmodul contributes a
// single oo-cusp term (m=2); multi-term bases (X0^21(2), the X0^10(N) pipeline) are still wrong and are
// gated off by the even-N guard in SchoferFormula. Do NOT extend this test to those bases until the
// vector-valued weight-3/2 Eisenstein coefficient is corrected at primes p | gcd(m, det).

procedure test_M0Multiplier()
    printf "Testing principled m=0 multiplier (Kudla-Yang) on X0^15(2)...";
    R<q> := LaurentSeriesRing(Rationals());

    // X0^15(2) star hauptmodul, principal part of the q-expansion at oo (cusp oo <-> eta = 0):
    //   2 q^-10 - 2 q^-2 + 2 q^-1   (constant/positive part irrelevant to the multiplier).
    foo := 2*q^-10 - 2*q^-2 + 2*q^-1;
    // The cusp-0 principal part distributes only to isotropic buckets whose local Whittaker vanishes
    // here, so it contributes 0; a trivial series suffices to pin the (oo-cusp) value. (The full cusp-0
    // path is covered end-to-end by SchoferIsometry.m.)
    f0 := R!1;

    Ldata := ShimuraCurveLattice(15, 2);
    mult := M0Multiplier(foo, f0, 15, 2, Ldata);

    // Validated value: sum_eta c_eta(0) = 4 (gives the +4 log 2 that matches all 19 Guo-Yang Table 45
    // discs, e.g. d = -7, -15, -60). It must be exactly the integer 4.
    assert mult eq 4;
    assert mult in Integers();

    printf "Done!\n";
end procedure;

test_M0Multiplier();
