// X_0^21(4)^* -> X_0^21(1)^*: degree 3, simple critical values at the star points of
// discriminants -4, -28, -84, -84' (branchdata.m), unramified at -7.  In the star coordinate x of
// models_21_1.m: x(-7) = 0, x(-28) = 1, x(-84) = (-5 +- sqrt(-3))/32, x(-4) = oo (the pole of both
// hauptmoduls).  Work in y = 1/x so that every critical value is finite: y(-4) = 0, y(-28) = 1,
// y(-84) = 8(-5 -+ sqrt(-3))/7, and the unramified -7 point sits at y = oo.
function degree3_maps(vals : F := Rationals())
    A<a1, a0, b1, b0, K, z> := PolynomialRing(F, 6);
    Pt<t> := PolynomialRing(A); Pv<v> := PolynomialRing(A);
    discs := [];
    for vv in [0, 1, -1, 2, -2] do
        cub := t^3 - vv*t^2 + (a1 - vv*b1)*t + (a0 - vv*b0);
        Append(~discs, Discriminant(cub));
    end for;
    pts := [0, 1, -1, 2, -2];
    Dv := Pv!0;
    for i in [1..5] do
        L := Pv!1;
        for j in [1..5] do
            if j ne i then L *:= (v - pts[j]) / (pts[i] - pts[j]); end if;
        end for;
        Dv +:= discs[i] * L;
    end for;
    target := K * &*[v - A!x : x in vals];
    eqs := [Coefficient(Dv, k) - Coefficient(target, k) : k in [0..4]] cat [K*z - 1];
    return ideal<A | eqs>;
end function;

Kf<s3> := QuadraticField(-3);
vals := [Kf!0, Kf!1, 8*(-5 - s3)/7, 8*(-5 + s3)/7];
I := degree3_maps(vals : F := Kf);
printf "dimension %o\n", Dimension(I);
V := Variety(I);
printf "solutions over Q(sqrt(-3)): %o\n", #V;
for p in V do
    printf "  R(t) = (t^3 + (%o) t + (%o)) / (t^2 + (%o) t + (%o))\n", p[1], p[2], p[3], p[4];
end for;
// the critical values of each solution, as a check, and the induced degree-3 map's ramification
Pz<T> := PolynomialRing(Kf);
for p in V do
    num := T^3 + p[1]*T + p[2]; den := T^2 + p[3]*T + p[4];
    cr := num*Derivative(den) - Derivative(num)*den;
    cvals := [];
    for rt in Roots(cr) do Append(~cvals, Evaluate(num, rt[1])/Evaluate(den, rt[1])); end for;
    printf "  critical values: %o   (wanted %o)\n", Sort(cvals), Sort(vals);
    printf "  simple critical points: %o of 4\n", #[rt : rt in Roots(cr) | rt[2] eq 1];
end for;
