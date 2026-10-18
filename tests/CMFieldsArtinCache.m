// tests/CMFieldsArtinCache.m
//
// Two discriminants with the same imaginary quadratic field K, computed one after the other while
// K is kept alive.  When the field of definition is K itself (the Atkin-Lehner group acts as all
// of Gal(H_R/K), or h(R) = 1 in the slow routine), the abelian extension's number field is K, and
// Magma's ArtinMap caches its result on that field: the second call would reuse the map of the
// first modulus and reject a prime dividing it ("Element is not in the codomain").  The failure
// appears only when K stays alive between the calls, which is why a local run can miss it.
//
// The expected fields are what the routines return when the two discriminants are computed in the
// opposite order in a fresh process (a regression value, not an independent verification); each
// is also checked against DegreeOfFieldOfDefinitionOfCMPoint.

procedure test_CMFieldsArtinCache()
    printf "Testing fields of definition at two discriminants sharing K, with K kept alive...";
    keep15 := NumberField(QuadraticOrder(BinaryQuadraticForms(-15)));
    X := CreateShimuraQuot(1, 30, {1, 2, 5, 10});
    for d in [-60, -15] do
        F := FieldOfDefinitionOfCMPoint(X, d);
        assert Type(F) eq FldRat and DegreeOfFieldOfDefinitionOfCMPoint(X, d) eq 1;
    end for;
    keep3 := NumberField(QuadraticOrder(BinaryQuadraticForms(-3)));
    Y := CreateShimuraQuot(10, 3, {1, 3, 10, 30});
    for d in [-27, -147] do
        F := FieldOfDefinitionOfCMPoint(Y, d);
        assert Type(F) eq FldRat and DegreeOfFieldOfDefinitionOfCMPoint(Y, d) eq 1;
    end for;
    // the slow routine, where h(R) = 1 puts the ring class field equal to K
    Z := CreateShimuraQuot(1, 39, {1, 3});
    for d in [-27, -3] do
        F := FieldOfDefinitionOfCMPointSlow(Z, d);
        assert AbsoluteDegree(F) eq 2 and IsIsomorphic(AbsoluteField(F), keep3);
        assert DegreeOfFieldOfDefinitionOfCMPoint(Z, d) eq 2;
    end for;
    assert keep15 eq keep15 and keep3 eq keep3;   // the fields stay alive to the end
    printf "Done!\n";
end procedure;
test_CMFieldsArtinCache();
