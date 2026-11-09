// Shape T of the certificate stage on X_0^6(35)/<w_3, w_5> (data/certificates/6_35.m): accepted, and refused, each
// for the reason it names, on the other quotient X_0^6(35)/<w_5, w_7> over the same E with the same branch points
// (a twist of it by rational 2-torsion, whose dual graph at 3 agrees with the square-root cover), over another
// elliptic quotient of X_0^6(35)/<w_3, w_5> where the dual graph at 3 does not separate the covers, and at p = 2,
// where the comparison does not apply.  A few seconds.

printf "Testing shape T of the certificate stage on X_0^6(35)/<w_3, w_5>...";
mk := function(W)
    Y := CreateShimuraQuot(6, 35, W); Y`g := GenusShimuraCurveQuotient(6, 35, W); Y`CurveID := 1;
    return Y;
end function;
Y := mk({Integers() | 1, 3, 5, 15});
assert Y`g eq 3;
cert := [c : c in ReadHyperellipticCertificates(6, 35) | c["shape"] eq "T"][1];
ok, why := VerifyHyperellipticCertificate(Y, cert);
assert ok and "HyperellipticCertificate T: not hyperelliptic" in why and "E(sqrt(c x)) has two singular points above it" in why;
Z := mk({Integers() | 1, 5, 7, 35});
bad := cert; bad["W"] := [1, 5, 7, 35]; bad["over"] := 3;
ok, why := VerifyHyperellipticCertificate(Z, bad);
assert not ok and "agrees with E(sqrt(c x))" in why;
bad := cert; bad["over"] := 2;
ok, why := VerifyHyperellipticCertificate(Y, bad);
assert not ok and "agrees with E(sqrt(c x))" in why;
bad := cert; bad["p"] := 2;
ok, why := VerifyHyperellipticCertificate(Y, bad);
assert not ok and "odd prime" in why;
printf " ok (accepted; refused for the twisted quotient over the same E, over another quotient, and at p = 2)\n";
