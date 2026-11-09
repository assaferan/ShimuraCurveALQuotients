// Shape D of the certificate stage on X_0^6(187)/<w_3, w_11, w_17> (data/certificates/6_187.m): accepted, and
// refused when E is another curve of the isogeny class 1122i of the star (its reduction types I_n are not the
// orders of the component groups); and, directly, the certified Heegner heights on that other curve are not
// lattice points, so no configuration survives there.  Several minutes.

printf "Testing shape D of the certificate stage on X_0^6(187)/<w_3, w_11, w_17>...";
H := {Integers() | 1, 3, 11, 17, 33, 51, 187, 561};
Y := CreateShimuraQuot(6, 187, H); Y`g := GenusShimuraCurveQuotient(6, 187, H); Y`CurveID := 1;
assert Y`g eq 3;
cert := ReadHyperellipticCertificates(6, 187)[1];
ok, why := VerifyHyperellipticCertificate(Y, cert);
assert ok and "HyperellipticCertificate D" in why;
others := [E2 : E2 in IsogenousCurves(EllipticCurve("1122i1")) | not IsIsomorphic(E2, EllipticCurve("1122i1"))];
assert #others ge 1;
bad := cert; bad["E"] := CremonaReference(others[1]);
ok, why := VerifyHyperellipticCertificate(Y, bad);
assert not ok and "reduction types" in why;
full := {Integers() | m : m in Divisors(6*187) | GCD(m, 6*187 div m) eq 1};
assert #ModelByHeegnerHeights(Y, full, MinimalModel(others[1]) : Certified := true) eq 0;
printf " ok (accepted; refused on %o, where no configuration survives)\n", CremonaReference(others[1]);
