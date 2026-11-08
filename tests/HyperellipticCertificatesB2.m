// Shape B2 of the certificate stage on X_0^6(341)/<w_3, w_11, w_62> (data/certificates/6_341.m): accepted,
// and refused, each for the reason it names, when the second curve is left out (the class has two
// principally polarised surfaces with real multiplication, so one curve does not decide), when it is
// replaced by its quadratic twist by -1 (its Hecke polynomials are not those of the star), and when its
// cited endomorphism is replaced by that of the first curve (it does not preserve the period lattice).
// About three minutes.

printf "Testing shape B2 of the certificate stage on X_0^6(341)/<w_3, w_11, w_62>...";
H := {Integers() | 1, 3, 11, 33, 62, 186, 682, 2046};
Y := CreateShimuraQuot(6, 341, H); Y`g := GenusShimuraCurveQuotient(6, 341, H); Y`CurveID := 1;
assert Y`g eq 4;
cert := ReadHyperellipticCertificates(6, 341)[1];
ok, why := VerifyHyperellipticCertificate(Y, cert);
assert ok and "HyperellipticCertificate B2" in why;
bad := cert; bad["stars"] := [cert["stars"][1]]; bad["rm"] := [cert["rm"][1]]; bad["conductor2"] := [2];
ok, why := VerifyHyperellipticCertificate(Y, bad);
assert not ok and "two non-isomorphic curves" in why;
bad := cert; bad["stars"][2] := <[-c : c in cert["stars"][2][1]], []>;   // y^2 = -f(x), the twist by -1
ok, why := VerifyHyperellipticCertificate(Y, bad);
assert not ok and "do not single out" in why;
bad := cert; bad["rm"][2] := cert["rm"][1];
ok, why := VerifyHyperellipticCertificate(Y, bad);
assert not ok and "period lattice" in why;
printf " ok (accepted; refused with one curve, with the twist by -1, and with the wrong endomorphism)\n";
