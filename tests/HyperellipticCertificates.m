// tests/HyperellipticCertificates.m
//
// The stage FilterByHyperellipticCertificates decides X_0^46(31)/W_full from its certificate
// data/certificates/46_31.m (shape A: four degree-two CM divisors related by Borcherds forms
// against a torsion bound of two; the inputs of tests/Hyperelliptic46_31.m), and refuses a
// certificate whose hypotheses fail, each control for the reason it names:
//   * the torsion bound from a single prime is far larger than 4, so the count fails;
//   * Z(-292) taken twice has degree four;
//   * a certificate whose group is not the curve's is for another curve.

printf "Testing the certificate stage on X_0^46(31)*...";
D := 46; N := 31;
full := {Integers() | m : m in Divisors(D*N) | GCD(m, (D*N) div m) eq 1};
X := CreateShimuraQuot(D, N, full);
X`g := GenusShimuraCurveQuotient(D, N, full); X`CurveID := 1;
assert X`g eq 3;

certs := ReadHyperellipticCertificates(D, N);
assert #certs eq 1;
cert := certs[1];

// controls first: they fail before the forms are rebuilt, except the count
bad := cert; bad["divisors"] := [cert["divisors"][i] : i in [1..3]] cat [[<-292, 2>]];
ok, why := VerifyHyperellipticCertificate(X, bad);
assert not ok and "degree 4" in why;
bad := cert; bad["W"] := [1, 2, 23, 46];
ok, why := VerifyHyperellipticCertificate(X, bad);
assert not ok and "another curve" in why;
star := rep{c : c in GetQuotientsAndGenera([X]) | IsStarCurve(c)};
assert JacobianTorsionBound(star : MaxPrime := 5, NumPrimes := 1) ge 4;   // else this control proves nothing
bad := cert; bad["torsion"] := [5, 1];
ok, why := VerifyHyperellipticCertificate(X, bad);
assert not ok and "divisors against a torsion bound" in why;

// the certificate itself, through the stage
curves := [X];
FilterByHyperellipticCertificates(~curves);
assert assigned curves[1]`IsSubhyp and curves[1]`IsSubhyp and curves[1]`IsHyp;
assert "HyperellipticCertificate A: 4 degree-two CM divisors related by 6 Borcherds forms, torsion bound 2" in curves[1]`TestInWhichProved;
// a decided curve is left alone
curves[1]`TestInWhichProved := "planted";
FilterByHyperellipticCertificates(~curves);
assert curves[1]`TestInWhichProved eq "planted";
printf " ok (accepted with torsion bound 2; refused for degree, for another curve, and for the count)\n";
