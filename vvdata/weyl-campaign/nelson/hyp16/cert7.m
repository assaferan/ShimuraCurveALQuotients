// The model predicts the two pairs of branch points are defined over Q(sqrt(13)) and Q(sqrt(-11)).
// The repository computes the field of definition of a CM point on a quotient directly.
AttachSpec("/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/ShimuraQuotients.spec");
SetColumns(0);
X := CreateShimuraQuot([* <"D",39>, <"N",11>, <"W", {IntegerRing()|1,11,39,429}>, <"g",3>, <"CurveID",0> *]);
printf "X = X_0^39(11)/{1,11,39,429}, genus %o\n\n", X`g;
for d in [-52, -132] do
    deg := DegreeOfFieldOfDefinitionOfCMPoint(X, d);
    printf "CM discriminant %o : degree of field of definition on the quotient = %o\n", d, deg;
    F := FieldOfDefinitionOfCMPoint(X, d);
    printf "   field = %o\n", F;
    if Type(F) eq FldNum or Type(F) eq FldQuad then
        printf "   degree %o, discriminant %o\n", Degree(F), Discriminant(F);
        if Degree(F) eq 2 then
            printf "   as Q(sqrt(m)) with m = %o\n", Squarefree(Abs(Discriminant(F)))*Sign(Discriminant(F));
        end if;
    end if;
end for;
printf "\nmodel predicts: disc -52 pair over Q(sqrt(13)), disc -132 pair over Q(sqrt(-11))\n";
print "CERT7 DONE";
