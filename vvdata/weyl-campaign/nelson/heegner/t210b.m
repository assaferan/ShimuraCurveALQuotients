AttachSpec("ShimuraQuotients.spec"); SetColumns(0);
star := func<a, b | a*b div GCD(a, b)^2>;
cases := [<1326, 1, {1,17,78,1326}>, <210, 13, {1,2,21,42,65,130,1365,2730}>];
for c in cases do
    D := c[1]; N := c[2]; W := {Integers() | w : w in c[3]};
    X := CreateShimuraQuot(D, N, W); X`g := GenusShimuraCurveQuotient(D, N, W);
    printf "\n=== X_0^%o(%o)/%o, genus %o ===\n", D, N, W, X`g;
    als := [d : d in Divisors(D*N) | GCD(d, (D*N) div d) eq 1 and d notin W];
    seen := {};
    for u in als do
        Wu := W join {star(u, w) : w in W};
        if Wu in seen then continue; end if; Include(~seen, Wu);
        if GenusShimuraCurveQuotient(D, N, Wu) ne 1 then continue; end if;
        printf "elliptic quotient by w_%o:\n", u;
        for v in Wu diff W do
            fp := NumFixedPointsByCMOrder(D, N, v);
            for d in Keys(fp) do
                if fp[d] eq 0 then continue; end if;
                K := FieldOfDefinitionOfCMPoint(X, d);
                dk := Type(K) eq FldRat select 1 else Discriminant(Integers(K));
                printf "   w_%o: %o fixed points, disc %o, field disc %o (degree %o)\n", v, fp[d], d, dk, Degree(K);
            end for;
        end for;
    end for;
end for;
