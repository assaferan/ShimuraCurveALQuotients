AttachSpec("ShimuraQuotients.spec");
SetColumns(0);
D := 15; N := 2;
Xstar := CreateShimuraQuot(D, N, Set(Divisors(D*N)));
Xstar`g := GenusShimuraCurveQuotient(D, N, Xstar`W); Xstar`CurveID := 0;
curves := GetQuotientsAndGenera([Xstar]);
_ := exists(star){c : c in curves | IsStarCurve(c)};
P<X> := PolynomialRing(Rationals());
Hs := AssociativeArray();
Hs[-588] := [X^3 - 191/54*X^2 + 343/432*X - 6889/48384, X^3 - 21403/4200*X^2 + 295733/75600*X - 6889/48384];
Hs[-1960] := [X^3 - 29036/6889*X^2 + 60676/6889*X - 40320/6889, X^3 + 1424209/41334*X^2 - 1416397/20667*X - 40320/6889];
for d in [-588, -1960] do
    printf "\n=== d = %o\n", d;
    flds := FieldsOfDefinitionOfCMPoint(star, d);
    for F in flds do
        if Type(F) eq FldRat then printf "  pipeline field: Q\n"; continue; end if;
        printf "  pipeline field: degree %o, discriminant %o, polredabs %o\n", Degree(F), Discriminant(MaximalOrder(F)), DefiningPolynomial(OptimizedRepresentation(F));
    end for;
    for H in Hs[d] do
        if not IsIrreducible(H) then printf "  H = %o REDUCIBLE: %o\n", H, Factorization(H); continue; end if;
        K := NumberField(H);
        printf "  H = %o: field discriminant %o, polredabs %o\n", H, Discriminant(MaximalOrder(K)), DefiningPolynomial(OptimizedRepresentation(K));
    end for;
end for;
