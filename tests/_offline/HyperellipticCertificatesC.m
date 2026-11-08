// Shape C of the certificate stage on X_0^15(38)/<w_2, w_15, w_19> (data/certificates/15_38.m): accepted, and
// refused, each for the reason it names, when the halves are another pairing of the four branch points (the
// dual test does not certify it), when E is another curve of the isogeny class (570e1, whose Tamagawa numbers
// are not the orders of the component groups), and when the twist argument is missing (E(Q)[2] = Z/2).  Needs
// the saved basis basis_1140_4000.txt in DualTestBasisDirectory(); about five minutes.

printf "Testing shape C of the certificate stage on X_0^15(38)/<w_2, w_15, w_19>...";
H := {Integers() | 1, 2, 15, 19, 30, 38, 285, 570};
Y := CreateShimuraQuot(15, 38, H); Y`g := GenusShimuraCurveQuotient(15, 38, H); Y`CurveID := 1;
assert Y`g eq 3;
cert := ReadHyperellipticCertificates(15, 38)[1];
ok, why := VerifyHyperellipticCertificate(Y, cert);
assert ok and "HyperellipticCertificate C" in why and "dual graph of Y at 3" in why;
bad := cert; bad["divisors"] := [ [<-228, 1>, <-12, 1>], [<-760, 1>, <-40, 1>] ];
ok, why := VerifyHyperellipticCertificate(Y, bad);
assert not ok and "does not certify" in why;
bad := cert; bad["E"] := "570e1";
ok, why := VerifyHyperellipticCertificate(Y, bad);
assert not ok and "Tamagawa" in why;
bad := AssociativeArray(); for k in Keys(cert) do if k ne "dualgraph" then bad[k] := cert[k]; end if; end for;
ok, why := VerifyHyperellipticCertificate(Y, bad);
assert not ok and "no twist argument" in why;
printf " ok (accepted; refused for another pairing, another curve of the class, and no twist argument)\n";
