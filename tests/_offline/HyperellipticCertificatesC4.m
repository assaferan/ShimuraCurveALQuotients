// Shape C4 of the certificate stage on X_0^6(323)/<w_6, w_19, w_34> (data/certificates/6_323.m): accepted, and
// refused, each for the reason it names, for a pair of discriminants that is not a genus character of -7752,
// and for a prime at which the Prym splits (5).  No basis needed; about ten minutes (identifying E by its
// traces up to the Sturm bound dominates).

printf "Testing shape C4 of the certificate stage on X_0^6(323)/<w_6, w_19, w_34>...";
H := {Integers() | 1, 6, 19, 34, 51, 114, 646, 969};
Y := CreateShimuraQuot(6, 323, H); Y`g := GenusShimuraCurveQuotient(6, 323, H); Y`CurveID := 1;
cert := ReadHyperellipticCertificates(6, 323)[1];
ok, why := VerifyHyperellipticCertificate(Y, cert);
assert ok and "C4" in why and "twisted modular symbol" in why;
bad := cert; bad["chi"] := <-3, 2584>;
ok, why := VerifyHyperellipticCertificate(Y, bad);
assert not ok and "genus character" in why;
bad := cert; bad["prym"] := 5;
ok, why := VerifyHyperellipticCertificate(Y, bad);
assert not ok and "reducible" in why;
printf " ok (accepted; refused for a non-genus character and for a prime where the Prym splits)\n";
