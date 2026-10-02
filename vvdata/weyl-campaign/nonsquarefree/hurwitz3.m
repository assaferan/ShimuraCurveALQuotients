// Degree-3 maps P^1 -> P^1 with four prescribed simple critical values.  Normal form (the Mobius
// freedom on the source used to put one preimage of oo at oo, kill t^2 and make the cubic monic):
//     R(t) = (t^3 + a1 t + a0) / (t^2 + b1 t + b0).
// v is a critical value iff disc_t(t^3 + a1 t + a0 - v (t^2 + b1 t + b0)) = 0, a quartic in v; it
// must equal K prod (v - v_i).  Solving the five coefficient equations gives every such R.
function degree3_maps(vals : F := Rationals())
    A<a1, a0, b1, b0, K, z> := PolynomialRing(F, 6);
    Pt<t> := PolynomialRing(A); Pv<v> := PolynomialRing(A);
    // discriminant of the cubic in t with v symbolic: compute as a polynomial in v by interpolation
    // at 5 values of v (the discriminant has degree 4 in v)
    discs := [];
    for vv in [0, 1, -1, 2, -2] do
        cub := t^3 - vv*t^2 + (a1 - vv*b1)*t + (a0 - vv*b0);
        Append(~discs, Discriminant(cub));
    end for;
    // Lagrange interpolation in v over A
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
    I := ideal<A | eqs>;
    return I;
end function;

// synthetic test: R = (t^3 + 2t - 1)/(t^2 + t + 3): recover it from its critical values
Q := Rationals(); Pq<x> := PolynomialRing(Q);
Pnum := x^3 + 2*x - 1; Qden := x^2 + x + 3;
crit := Pnum*Derivative(Qden) - Derivative(Pnum)*Qden;   // negative of the usual, same roots
Kc<r> := SplittingField(crit);
Pk := PolynomialRing(Kc);
cvals := [Evaluate(Pk!Pnum, rt[1]) / Evaluate(Pk!Qden, rt[1]) : rt in Roots(Pk!crit)];
printf "synthetic critical values: %o\n", cvals;
I := degree3_maps(cvals : F := Kc);
printf "dimension of the solution set: %o\n", Dimension(I);
V := Variety(I);
printf "number of degree-3 maps with these critical values: %o\n", #V;
printf "the synthetic map recovered: %o\n", exists{p : p in V | p[1] eq 2 and p[2] eq -1 and p[3] eq 1 and p[4] eq 3};
