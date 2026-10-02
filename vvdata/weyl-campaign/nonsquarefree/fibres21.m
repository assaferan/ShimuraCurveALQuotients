// The degree-3 map X_0^21(4)^* -> X_0^21(1)^* found by the Hurwitz solve, and its fibres over the
// four CM points of the base.  Coordinate on the base: y = 1/s, with s the hauptmodul of
// models_21_1.m (zero at d = -7, pole at -4, and s + s~ = 1 gives s(-28) = 1); so
// y(-4) = 0, y(-28) = 1, y(-7) = oo, y(-84) = 8(-5 -+ sqrt(-3))/7.
Q := Rationals(); P<t> := PolynomialRing(Q);
num := t^3 - 4/3*t + 16/27; den := t^2 + 29/12*t + 22/9;
printf "R = (%o) / (%o)\n", num, den;
Kf<s3> := QuadraticField(-3);
named := [<-4, Kf!0>, <-28, Kf!1>, <-84, 8*(-5 - s3)/7>, <-84, 8*(-5 + s3)/7>];
PK<T> := PolynomialRing(Kf);
for pr in named do
    f := PK!num - pr[2]*PK!den;
    printf "fibre over d = %-4o (y = %o): %o\n", pr[1], pr[2], [<Degree(fa[1]), fa[2]> : fa in Factorization(f)];
end for;
printf "fibre over d = -7 (y = oo): t = oo, and %o\n", [<Degree(fa[1]), fa[2]> : fa in Factorization(den)];
printf "  the two finite ones generate %o (want Q(sqrt(-7)))\n", Discriminant(NumberField(P!den));
// the critical points, exactly
cr := num*Derivative(den) - Derivative(num)*den;
printf "critical points: %o\n", [<Degree(fa[1]), fa[2], Coefficients(fa[1])> : fa in Factorization(cr)];
for fa in Factorization(cr) do
    if Degree(fa[1]) eq 1 then
        r := -Coefficient(fa[1],0)/Coefficient(fa[1],1);
        printf "  t = %o  ->  y = %o\n", r, Evaluate(num, r)/Evaluate(den, r);
    else
        Kl := NumberField(fa[1]); printf "  a conjugate pair over disc %o\n", Discriminant(Kl);
    end if;
end for;
