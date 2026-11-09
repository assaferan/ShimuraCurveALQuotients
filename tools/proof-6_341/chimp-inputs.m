// The two inputs of data/certificates/6_341.m that the certificate stage cites rather than recomputes,
// reproduced with CHIMP: for both curves a certified endomorphism over Q whose tangent matrix is the
// stored one, and the proved conductor exponent at 2.  About a minute.
//   magma -b CHIMP:=/path/to/CHIMP.spec tools/proof-6_341/chimp-inputs.m < /dev/null   (from the repo root)
AttachSpec("ShimuraQuotients.spec");
AttachSpec(CHIMP);
cert := ReadHyperellipticCertificates(6, 341)[1];
PQ<x> := PolynomialRing(Rationals());
for i -> c in cert["stars"] do
    C := HyperellipticCurve(PQ!c[1], PQ!c[2]);
    f := HyperellipticPolynomials(SimplifiedModel(C));
    F := RationalsExtra(100);
    _, certs := CertifiedEndomorphismAlgebra(HyperellipticCurve(ChangeRing(f, F)));
    T := [Rationals() | a : a in cert["rm"][i]];
    ok1 := exists{e : e in certs | [Rationals() | a : a in Eltseq(e[1])] eq T};
    e2 := EvenConductorExponent_Genus2(C, 2 : Proof := true);
    printf "curve %o: certified endomorphism with tangent matrix %o: %o; proved conductor exponent at 2: %o (stored %o)\n",
        i, T, ok1, e2, cert["conductor2"][i];
    assert ok1 and e2 eq cert["conductor2"][i];
end for;
print "CHIMP inputs reproduced";
quit;
