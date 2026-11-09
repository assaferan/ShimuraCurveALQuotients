// Shape O of the certificate stage on X_0^158(3)/<w_2, w_237> (data/certificates/158_3.m): accepted, and
// refused, each for the reason it names, when E is a curve of another isogeny class of conductor 158, and when
// p is 2, which divides D.  The mechanism is checked on the height of the chi-component, chi the genus
// character of Q(sqrt -79) on discriminant -632: on Z = X_0^158(1)/<w_2>, which the dual graph at 2 shows is
// not hyperelliptic (recorded in data/curves_after_UpdateCurves8.dat), Cai-Shu-Tian's formula gives a nonzero
// height, and with the form old at 3 the factor (1 + chi(p_3))^2 makes it 0.  A few seconds.

printf "Testing shape O of the certificate stage on X_0^158(3)/<w_2, w_237>...";
H := {Integers() | 1, 2, 237, 474};
Y := CreateShimuraQuot(158, 3, H); Y`g := GenusShimuraCurveQuotient(158, 3, H); Y`CurveID := 1;
assert Y`g eq 3;
cert := ReadHyperellipticCertificates(158, 3)[1];
ok, why := VerifyHyperellipticCertificate(Y, cert);
assert ok and "HyperellipticCertificate O" in why;
bad := cert; bad["E"] := "158b1";
ok, why := VerifyHyperellipticCertificate(Y, bad);
assert not ok and "differ at" in why;
bad := cert; bad["p"] := 2;
ok, why := VerifyHyperellipticCertificate(Y, bad);
assert not ok and "prime to D" in why;
E := EllipticCurve("158a1");
cZ := HeegnerHeightComponents(E, 158, 1, 1, -632 : Certified := true);
cX := HeegnerHeightComponents(E, 158, 3, 1, -632 : Certified := true, OldformPrime := 3);
hZ := [x[6] : x in cZ | <x[1], x[2]> eq <-79, 8>][1];
hX := [x[6] : x in cX | <x[1], x[2]> eq <-79, 8>][1];
assert hZ gt 1 and hX eq 0;
printf " ok (accepted; refused for another class and for p dividing D; chi-height %o on Z, 0 old at 3)\n", RealField(6)!hZ;
