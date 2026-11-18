// tests/HyperellipticCertificates.m
//
// The stage FilterByHyperellipticCertificates decides X_0^46(31)/W_full from its certificate
// data/certificates/46_31.m (shape A: four degree-two CM divisors related by Borcherds forms
// against a torsion bound of two; the inputs of tests/Hyperelliptic46_31.m), and refuses a
// certificate whose hypotheses fail, each control for the reason it names:
//   * the torsion bound from a single prime is far larger than 4, so the count fails;
//   * Z(-292) taken twice has degree four;
//   * a certificate whose group is not the curve's is for another curve.
// Shape B on X_0^2030(1)/<w_10, w_14, w_58> (genus four over the genus-two star, branched at
// P_-8 + P_-280; data/certificates/2030_1.m): accepted, and refused when the first divisor is not the
// branch divisor, and when the torsion bound, from a single prime, is not 1.  The acceptance names
// p = 3 as a Frobenius of order 5 on J[2] of the star: T_3 on the newforms of level 2030 with every
// w_q = -1 has characteristic polynomial x^2 + 3x + 1 (Magma modular symbols), so a = 3, b = 2*3 + 1.
// Shape N on X_0^141(2)/{1, 6, 94, 141} (data/certificates/141_2.m): proved not hyperelliptic, the stage
// records it so, and the certificate is refused when its only prime is 47, where the class is zero.

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

// shape B
H := {Integers() | 1, 10, 14, 35, 58, 145, 203, 2030};
Y := CreateShimuraQuot(2030, 1, H); Y`g := GenusShimuraCurveQuotient(2030, 1, H); Y`CurveID := 2;
assert Y`g eq 4;
certB := ReadHyperellipticCertificates(2030, 1)[1];
ok, why := VerifyHyperellipticCertificate(Y, certB);
assert ok and "HyperellipticCertificate B" in why and "related to 2 other divisor(s)" in why;
assert "Frobenius of order 5 on J[2] at p = 3" in why;
bad := certB; bad["divisors"] := [certB["divisors"][2], certB["divisors"][1], certB["divisors"][3]];
ok, why := VerifyHyperellipticCertificate(Y, bad);
assert not ok and "branch divisor" in why;
starB := rep{c : c in GetQuotientsAndGenera([Y]) | IsStarCurve(c)};
assert JacobianTorsionBound(starB : MaxPrime := 5, NumPrimes := 1, Check := false) gt 1;   // else this control proves nothing
bad := certB; bad["torsion"] := [5, 1];
ok, why := VerifyHyperellipticCertificate(Y, bad);
assert not ok and "not 1" in why;
// shape N
HN := {Integers() | 1, 6, 94, 141};
Z := CreateShimuraQuot(141, 2, HN); Z`g := GenusShimuraCurveQuotient(141, 2, HN); Z`CurveID := 3;
certN := [c : c in ReadHyperellipticCertificates(141, 2) | c["shape"] eq "N"][1];
ok, why := VerifyHyperellipticCertificate(Z, certN);
assert ok and "not hyperelliptic" in why;
zs := [Z]; FilterByHyperellipticCertificates(~zs);
assert assigned zs[1]`IsSubhyp and not zs[1]`IsSubhyp and not zs[1]`IsHyp;
bad := certN; bad["phi"] := [47];
ok, why := VerifyHyperellipticCertificate(Z, bad);
assert not ok and "class zero" in why;
printf " ok (A: accepted with torsion bound 2; refused for degree, for another curve, and for the count; B: accepted, refused for the branch divisor and the torsion; N: accepted, recorded as not hyperelliptic, refused at 47)\n";
