// Shape V of the certificate stage (data/certificates/6_35.m): accepted on both genus-three quotients of X_0^6(35)
// over the same elliptic quotient E = X_0^6(35)/<w_3, w_5, w_7>, and refused, each for the reason it names, on two
// curves that the certificate stage proves hyperelliptic by shape D (X_0^210(13)/W and X_0^1326(1)/W of
// data/certificates/210_13.m, 1326_1.m), on a curve whose intermediate quotients are not all elliptic, and when the
// involution named is in W.  Each of the two hyperelliptic curves meets every other hypothesis over one of its
// elliptic quotients, where its two branch places have the same field, so that hypothesis cannot be dropped.  A few
// seconds.

printf "Testing shape V of the certificate stage...";
mk := function(D, N, W)
    Y := CreateShimuraQuot(D, N, {Integers() | m : m in W}); Y`g := GenusShimuraCurveQuotient(D, N, Y`W); Y`CurveID := 1;
    return Y;
end function;
vcert := function(D, N, W, m)
    c := AssociativeArray(); c["D"] := D; c["N"] := N; c["W"] := W; c["shape"] := "V"; c["over"] := m; c["file"] := "test";
    return c;
end function;
certs := [c : c in ReadHyperellipticCertificates(6, 35) | c["shape"] eq "V"];
assert #certs eq 2;
for c in certs do
    ok, why := VerifyHyperellipticCertificate(mk(6, 35, c["W"]), c);
    assert ok and "HyperellipticCertificate V: not hyperelliptic" in why and "[ -420, -84 ]" in why;
end for;
// hyperelliptic curves, refused over each of their three elliptic quotients
for t in [<210, 13, [1, 2, 21, 42, 65, 130, 1365, 2730]>, <1326, 1, [1, 17, 78, 1326]>] do
    D, N, W := Explode(t);
    Y := mk(D, N, W); L := D*N;
    same := false;
    for m in [m : m in Divisors(L) | GCD(m, L div m) eq 1 and m notin W] do
        ok, why := VerifyHyperellipticCertificate(Y, vcert(D, N, W, m));
        assert not ok and ("are not two places of degree two" in why or "the same field of definition" in why);
        same or:= "the same field of definition" in why;
    end for;
    assert same;
end for;
ok, why := VerifyHyperellipticCertificate(mk(6, 23, [1, 2]), vcert(6, 23, [1, 2], 3));
assert not ok and "not 1, 1, 1" in why;
ok, why := VerifyHyperellipticCertificate(mk(6, 35, [1, 3, 5, 15]), vcert(6, 35, [1, 3, 5, 15], 5));
assert not ok and "outside W" in why;
printf " ok (both quotients of X_0^6(35) accepted; refused on two hyperelliptic curves, one with equal fields over each, for a non-elliptic intermediate quotient, and inside W)\n";
