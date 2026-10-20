// The equation of a hyperelliptic star quotient from two of its functions, given by their values at
// the CM points and their divisors.  This is the last step of Guo-Yang's recipe (Compositio 153
// (2017), Examples 42 and 43): having found Borcherds forms with prescribed CM divisors and
// evaluated them by Schofer's formula, one fits the relation between the two coordinates they define
// and reads the curve off it.  At genus one the relation has bidegree (2, 2) and the curve is read
// by Magma's EllipticCurve (GenusOneStar.m); here x is a function of degree two, so the relation is
// quadratic in y and completing the square gives a model y^2 = f(x) at any genus.
//
// WHY THE VALUES ALONE DO NOT SUFFICE.  A relation of bidegree (n, 2) has 3(n + 1) coefficients, so
// pinning it up to a scalar takes 3n + 2 conditions, while a star carries only as many rational CM
// points as it carries -- ten at X_0^302(1)*, against the eleven that its bidegree (3, 2) relation
// would need.  What closes the gap is the divisor of y, which is known exactly because it was
// prescribed to the Borcherds search: writing the relation as A(x) y^2 + B(x) y + C(x) = 0, the
// quotients -B/A and C/A are the trace and the norm of y over Q(x), so the roots of A are the
// images under x of the poles of y and the roots of C are the images of its zeros.  The
// multiplicity of a root x_0 is the sum of the orders of y over the whole fibre of x above x_0,
// which is the same number whether or not the fibre is ramified, and the points of that fibre carry
// a nonzero coefficient of div y, so they are all listed in it.

function is_infinity(v)
    return Type(v) eq Infty;
end function;

// The relation of bidegree (dx, dy) whose coefficients, listed as x^a y^b with a running fastest
// over [0..dx] and b over [0..dy], satisfy the linear conditions in rows.
function relation_from_conditions(rows, dx, dy)
    K := Kernel(Transpose(Matrix(Rationals(), rows)));
    error if Dimension(K) ne 1,
        Sprintf("the conditions leave a %o-dimensional space of relations", Dimension(K));
    kv := Basis(K)[1];
    P<x, y> := PolynomialRing(Rationals(), 2);
    mons := [x^a*y^b : a in [0..dx], b in [0..dy]];
    F := &+[kv[i]*mons[i] : i in [1..#mons]];
    F *:= LCM([Denominator(c) : c in Coefficients(F)]);
    F /:= GCD([Integers() | c : c in Coefficients(F)]);
    if LeadingCoefficient(F) lt 0 then F := -F; end if;
    return F;
end function;

intrinsic BidegreeRelation(xs::List, ys::List, dx::RngIntElt, dy::RngIntElt) -> RngMPolElt
{The relation F(x, y) = 0 of bidegree (dx, dy) satisfied by the pairs (xs[i], ys[i]), where
Infinity() marks a pole, when the pairs determine it up to a scalar.  F is returned with integer
coefficients of content one and positive leading coefficient.}
    require #xs eq #ys : "xs and ys must have the same length";
    require dx ge 1 and dy ge 1 : "both degrees must be at least 1";
    rows := [];
    for i in [1..#xs] do
        X := is_infinity(xs[i]) select [Rationals() | 0, 1] else [Rationals() | 1, xs[i]];   // (X0 : X1), x = X1/X0
        Y := is_infinity(ys[i]) select [Rationals() | 0, 1] else [Rationals() | 1, ys[i]];
        Append(~rows, [X[1]^(dx-a)*X[2]^a*Y[1]^(dy-b)*Y[2]^b : a in [0..dx], b in [0..dy]]);
    end for;
    return relation_from_conditions(rows, dx, dy);
end intrinsic;

intrinsic HyperellipticRelation(ds::SeqEnum[RngIntElt], xs::List, ys::List, divy::SeqEnum) -> RngMPolElt
{The relation A(x) y^2 + B(x) y + C(x) = 0 between a function x of degree two on a star quotient and
a second function y, from their values xs, ys at the rational CM points of discriminants ds
(Infinity() for a pole) and the divisor divy of y, given as pairs <d, a_d> with d among ds.  Besides
the values at the points, the fit uses the roots that the divisor forces on A and on C: see the
header of this file.  Every point of divy must be rational, so that its coefficient is its order.}
    require #xs eq #ds and #ys eq #ds : "one value of each coordinate per discriminant";
    require forall{c : c in divy | c[1] in ds} : "every point of div y must be among ds";
    require &+[c[2] : c in divy] eq 0 : "div y must have degree zero";
    n := &+[c[2] : c in divy | c[2] gt 0];
    val := func<d | ys[Index(ds, d)]>;
    xval := func<d | xs[Index(ds, d)]>;
    for c in divy do
        require (c[2] gt 0 and val(c[1]) eq 0) or (c[2] lt 0 and is_infinity(val(c[1]))) :
            Sprintf("the value of y at %o is %o, which does not match its multiplicity %o in div y",
                    c[1], val(c[1]), c[2]);
    end for;

    // the roots the divisor forces, with multiplicity: those of A from the poles of y, those of C
    // from its zeros, each grouped over the fibre of x and with the points above x = infinity left out
    forced := function(part)
        rts := AssociativeArray();
        for c in [c : c in divy | part*c[2] gt 0] do
            x0 := xval(c[1]);
            if is_infinity(x0) then continue; end if;
            rts[x0] := (IsDefined(rts, x0) select rts[x0] else 0) + part*c[2];
        end for;
        return rts;
    end function;

    rows := [];
    for j->d in ds do
        X := is_infinity(xs[j]) select [Rationals() | 0, 1] else [Rationals() | 1, xs[j]];
        Y := is_infinity(ys[j]) select [Rationals() | 0, 1] else [Rationals() | 1, ys[j]];
        Append(~rows, [X[1]^(n-a)*X[2]^a*Y[1]^(2-b)*Y[2]^b : a in [0..n], b in [0..2]]);
    end for;
    // A is the coefficient of y^2 and vanishes where y has a pole, C that of y^0 and vanishes where
    // y has a zero; a root of multiplicity k at x_0 is the vanishing of the polynomial and of its
    // first k-1 derivatives there
    for pr in [<2, -1>, <0, 1>] do
        b := pr[1];
        rts := forced(pr[2]);
        for x0 in Keys(rts) do
            for e in [0..rts[x0] - 1] do
                row := [Rationals() | 0 : a in [0..n], bb in [0..2]];
                for a in [e..n] do
                    row[(a + 1) + b*(n + 1)] := Factorial(a) / Factorial(a - e) * x0^(a - e);
                end for;
                Append(~rows, row);
            end for;
        end for;
    end for;
    F := relation_from_conditions(rows, n, 2);
    require Degree(F, 2) eq 2 : "the relation found is not quadratic in y";
    return F;
end intrinsic;

intrinsic HyperellipticModelFromRelation(F::RngMPolElt) -> CrvHyp, RngMPolElt, RngMPolElt
{A hyperelliptic model Y^2 = f(x) of the curve F(x, y) = 0, which must be quadratic in y, together
with the numerator and denominator of Y as a function of x and y.  Writing F = A y^2 + B y + C, the
substitution Y = (2 A y + B)/s with s^2 the square part of the discriminant B^2 - 4 A C generates the
same field over Q(x), so the model is the same curve; its genus is that of F = 0 when F is
irreducible and the curve has no singularity off the plane.}
    P := Parent(F);
    require Rank(P) eq 2 : "F must be a polynomial in two variables";
    require Degree(F, 2) eq 2 : "F must be quadratic in its second variable";
    x := P.1; y := P.2;
    R<t> := PolynomialRing(Rationals());
    coeff := function(j)
        c := R!0;
        for m in Monomials(F) do
            if Degree(m, 2) eq j then c +:= MonomialCoefficient(F, m) * t^Degree(m, 1); end if;
        end for;
        return c;
    end function;
    Ap := coeff(2); Bp := coeff(1); Cp := coeff(0);
    disc := Bp^2 - 4*Ap*Cp;
    require disc ne 0 : "the relation is a square in y, so it defines no double cover";
    // divide out the square part, which only rescales Y; a square constant factor counts too
    s := R!1;
    for fa in Factorization(disc) do
        s *:= fa[1]^(fa[2] div 2);
    end for;
    f := ExactQuotient(disc, s^2);
    sq, r := IsSquare(LeadingCoefficient(f));
    if sq then f /:= LeadingCoefficient(f); s *:= r; end if;
    require Degree(f) ge 3 : Sprintf("the squarefree part of the discriminant has degree %o, too small for a hyperelliptic curve", Degree(f));
    H := HyperellipticCurve(f);
    num := Evaluate(2*Ap, x)*y + Evaluate(Bp, x);
    return H, num, Evaluate(s, x);
end intrinsic;
