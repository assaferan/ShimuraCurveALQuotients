// Preconditions of ModelByHeegnerHeights at 39_11 and 55_7: elliptic quotients by signs, branch
// discriminants and fields, class groups, ranks of E and of its twists by the branch fields.
AttachSpec("ShimuraQuotients.spec"); SetColumns(0);
for c in [<39, 11, {1,11,39,429}>] do
    D, N, W := Explode(c); W := {Integers() | w : w in W};
    X := CreateShimuraQuot(D, N, W); X`g := GenusShimuraCurveQuotient(D, N, W);
    printf "\n=== X_0^%o(%o)/%o, genus %o ===\n", D, N, W, X`g;
    for q in EllipticQuotientsBySigns(D, N, W) do
        Wp, E := Explode(q);
        printf "quotient by %o: class of %o (rank %o, torsion %o)\n", Setseq(Wp diff W), CremonaReference(E), Rank(E), #TorsionSubgroup(E);
        for b in BranchPointData(X, Wp) do
            F := b[3];
            if Type(F) eq FldRat then printf "   w_%o: disc %o, RATIONAL points (%o on X)\n", b[1], b[2], b[4]; continue; end if;
            dF := Discriminant(Integers(F)); Cl := ClassGroup(QuadraticField(b[2]));
            printf "   w_%o: disc %o (h = %o, exponent %o, fundamental %o), field disc %o, %o points on X; twist rank %o\n",
                b[1], b[2], #Cl, Exponent(Cl), IsFundamentalDiscriminant(b[2]), dF, b[4], Rank(QuadraticTwist(E, dF));
        end for;
    end for;
end for;
