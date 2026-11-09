// SCHOFER'S FORMULA ON A DUAL CERTIFICATE AT A NON-FUNDAMENTAL DISCRIMINANT.  About two minutes.
//
// data/dual/cert_2002_1.m is the principal part at the cusp oo of a weight-1/2 form f on Gamma_0(4004)
// holomorphic at every other cusp, whose Borcherds lift psi has divisor lambda (P_-11 + P_-91 - P_-148
// - P_-232) on X_0^2002(1)/<w_2,w_13,w_77>, lambda = 222706250475728656 (tests/DualCertificate2002_1_divisor.m).
// The genus-two star has a coordinate x with divisor P_-148 + P_-232 - P_-11 - P_-91 and x(-99) = 1/3,
// x(-280) = -1 (paper/star-2002-proof.tex, the proof of the model), so psi = c x^(-lambda) and
//     log|psi(tau_-99)| - log|psi(tau_-280)| = lambda log 3.
// The discriminant -99 = 9 * (-11) has conductor 3, prime to 2002.  Schofer's formula from the expansion
// at oo alone -- the sum over m >= 1 -- gives the value at -280 and at every fundamental discriminant
// but misses, at -99, the m = 0 term at the conductor prime: a fibre sum over the pole of f at -11,
// which is exactly lambda log 3 here (M0FibreCorrection with an empty cusp-0 side, since f is holomorphic
// at the cusp 0).  This test checks the difference, and that the fibre term alone is lambda log 3.

// the certificate is read at top level: evaluating a data file inside a procedure crashes Magma here
cert := eval (Read("data/dual/cert_2002_1.m") cat " return cert;");

procedure test_certificate_schofer_2002_1(cert)
    printf "Schofer's formula on the 2002 certificate: log|psi(-99)| - log|psi(-280)| = lambda log 3...";
    D := 2002; N := 1;
    lambda := 222706250475728656;
    R<q> := LaurentSeriesRing(Rationals());
    B := #cert;
    f := &+[ R | cert[n] * q^(-n) : n in [1..B] | cert[n] ne 0 ] + O(q^2);
    assert cert[11] eq -lambda;                      // the pole at -11 = -99/9 (index n = 11 = 4*11/2^2)
    Ldata := ShimuraCurveLattice(D, N);
    Q := Ldata`Q;
    lam99 := ElementOfNorm(Q, 99, Ldata`O, Ldata`basis_L);
    corr := M0FibreCorrection([f], [BigO(q)], [Rationals() | 0], -99, 3, Q : Lambda := lam99, M := 2*D*N, Unimodular := true);
    assert corr[1] eq lambda;
    assert ScaleForSchofer(-99, D, N) eq -1/4;       // one CM point of discriminant -99 on the star
    v99 := SchoferFormula(f, -99, D, N, Ldata : Lambda := lam99);
    v280 := SchoferFormula(f, -280, D, N, Ldata);
    assert v99 - v280 eq LogSum(Rationals()!lambda, 3);
    // without the fibre term the value at -99 is that of a fundamental point over x = 1 (wrong)
    bare := SchoferFormula(f, -99, Q, lam99, ScaleForSchofer(-99, D, N));
    assert v99 - bare eq LogSum(Rationals()!lambda, 3);
    printf " Done!\n";
end procedure;

test_certificate_schofer_2002_1(cert);
