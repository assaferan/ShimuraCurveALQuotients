// X_0^6(379)* IS HYPERELLIPTIC, PART TWO: THE FORM EXISTS.  About an hour.
//
// A weakly holomorphic form of weight 1/2 on Gamma_0(4548) with the theta multiplier, holomorphic
// away from the cusp infinity and with principal part sum_n c(-n) q^-n there, exists if and only if
// sum_n c(-n) a_g(n) = 0 for every cusp form g of weight 3/2 on Gamma_0(4548) with the theta^3
// multiplier (the residue theorem for the weight-two differential f g dtau on X_0(4548), and Serre
// duality for the converse).  This test computes that space -- its dimension, 377, is the
// obstruction dimension the Borcherds pipeline predicts at level 4548 -- to precision 700 and
// checks the pairing of the banked principal part with every basis form.  The divisor of the lift
// is tests/DualCertificate6_379_divisor.m.

// the certificate is read at top level: evaluating a data file inside a procedure crashes Magma here
cert := eval (Read("data/dual/cert_6_379.m") cat " return cert;");

procedure test_dual_certificate_6_379(cert)
    printf "a weight-1/2 form on Gamma_0(4548) with the banked principal part exists...";
    B := #cert;
    S := CuspidalSubspace(HalfIntegralWeightForms(Gamma0(4548), 3/2));
    assert Dimension(S) eq 377;
    Bas := Basis(S, B + 2);
    pair := [ &+[ cert[n] * Coefficient(g, n) : n in [1..B] ] : g in Bas ];
    assert forall{ x : x in pair | x eq 0 };
    // and a principal part that is NOT a certificate fails: flip the sign of one coefficient
    c2 := cert; n0 := Minimum([n : n in [1..B] | cert[n] ne 0]); c2[n0] := -c2[n0];
    pair2 := [ &+[ c2[n] * Coefficient(g, n) : n in [1..B] ] : g in Bas ];
    assert exists{ x : x in pair2 | x ne 0 };
    printf " [377 cusp forms, pairing zero; the perturbed part fails] Done!\n";
end procedure;

test_dual_certificate_6_379(cert);
