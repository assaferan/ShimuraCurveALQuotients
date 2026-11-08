// Shape C on the two quotients whose star's form is old at 3, X_0^14(51)/<w_6, w_14, w_34> and
// X_0^38(15)/<w_3, w_10, w_38> (data/certificates/14_51.m, 38_15.m): both accepted, and the first refused,
// each for the reason it names, for another pairing of its four rational branch points (the dual test does not
// certify it) and for E = 238a2, the other curve of the class (its reduction types I_n do not have the orders
// of the component groups).  Needs basis_1428_4000.txt and basis_1140_4000.txt in DualTestBasisDirectory();
// about four minutes.

printf "Testing shape C on X_0^14(51)/<w_6, w_14, w_34> and X_0^38(15)/<w_3, w_10, w_38>...";
H := {Integers() | 1, 6, 14, 21, 34, 51, 119, 714};
Y := CreateShimuraQuot(14, 51, H); Y`g := GenusShimuraCurveQuotient(14, 51, H); Y`CurveID := 1;
assert Y`g eq 3;
cert := ReadHyperellipticCertificates(14, 51)[1];
ok, why := VerifyHyperellipticCertificate(Y, cert);
assert ok and "HyperellipticCertificate C" in why;
bad := cert; bad["divisors"] := [ [<-1428, 1>, <-168, 1>], [<-408, 1>, <-8, 1>] ];
ok, why := VerifyHyperellipticCertificate(Y, bad);
assert not ok and "does not certify" in why;
bad := cert; bad["E"] := "238a2";
ok, why := VerifyHyperellipticCertificate(Y, bad);
assert not ok and "reduction types" in why;
H := {Integers() | 1, 3, 10, 30, 38, 95, 114, 285};
Y := CreateShimuraQuot(38, 15, H); Y`g := GenusShimuraCurveQuotient(38, 15, H); Y`CurveID := 1;
assert Y`g eq 3;
cert := ReadHyperellipticCertificates(38, 15)[1];
ok, why := VerifyHyperellipticCertificate(Y, cert);
assert ok and "HyperellipticCertificate C" in why;
printf " ok (both accepted; 14_51 refused for another pairing and for 238a2)\n";
