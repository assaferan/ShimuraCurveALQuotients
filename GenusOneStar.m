// A model of a star quotient of genus one, following Guo-Yang, "Equations of hyperelliptic
// Shimura curves", Compositio Math. 153 (2017), Example 42.
//
// When the star curve X has genus one no divisor of degree one is principal, so the hauptmodul
// search of BorcherdsForms cannot succeed; what exists instead are functions of degree two.  A
// function x with a double pole at a CM point P is pinned by two Borcherds forms psi_0, psi_1 with
// divisor (two CM points) - 2P, obtained from BorcherdsForms(... : Divisors := ...): x = c_0 psi_0
// vanishes at the zeros of psi_0, and 1 - x = c_1 psi_1 vanishes at those of psi_1.  Schofer's
// formula gives |psi_j| at every CM point, which fixes |c_0|, |c_1| and, through |x| and |1 - x|
// together, the sign of x.  Two such functions with distinct poles generate the function field, and
// the relation between them, of bidegree (2, 2), is read off their values at nine or more CM points.
// tests/GuoYangExample42.m runs the whole recipe on X_0^142(1)* against the paper's own numbers.
//
// Where Guo-Yang's Assumption 40 enters.  A Borcherds form carries a character of the group that
// uniformises the star quotient (Guo-Yang, Lemma 5 and Remark 41), and it is a function on the
// star only when that character is trivial, which Guo-Yang assume (their Assumption 40; Example
// 42 is stated under it).  The search here does not check the character: it only finds a form
// whose divisor on the star is the prescribed one.
// When the principal part of the input is INTEGRAL, that character has finite order (Borcherds,
// "The Gross-Kohnen-Zagier theorem in higher dimensions", Duke Math. J. 97 (1999), proof of
// Theorem 4.1, with the correction in Duke Math. J. 105 (2000), 183-184 = arXiv alg-geom/9710002v3:
// the printed 1999 version omits the finite-order condition from the definition of a principal
// Heegner divisor, and the correction restores it and fills the gap in the proof); the search is
// run with IntegralSolution so that it is.  Then some power of the
// form is a function, so the prescribed divisor, which is rational, has a class of finite order
// in E(Q); and the order of the character is exactly the order of that class, since psi^n is a
// function precisely when n times the divisor is principal.
// When E(Q) has no torsion the class is therefore trivial, a function g with that divisor exists,
// and the form divided by g has no zeros or poles, so the character is trivial and the assumption
// holds: this is the case at X_0^142(1)* (142a1 has trivial torsion).  When E(Q) has rational
// torsion the assumption is a genuine hypothesis -- the form may be a root of a function, and its
// Schofer values are then those of the root -- and a model built from such forms is conditional
// on it.  Replacing each form by its e-th power, e the exponent of the torsion of E(Q), gives
// functions for certain, with the e-th powers of the values.

import "BorcherdsForms.m" : basis_of_weakly_holomorphic_forms, coeffs_to_divisor_matrix;

function is_infinite(v)
    return Type(v) eq Infty;
end function;

intrinsic BorcherdsFormValues(fs::SeqEnum[EtaQuot], ds::SeqEnum[RngIntElt], Xstar::ShimuraQuot : Degrees := []) -> SeqEnum
{The absolute values |psi_f| of the Borcherds forms of the forms in fs at the CM points of Xstar of
discriminants ds, one list per form (a zero gives 0, a pole Infinity()).  Degrees gives the degree
of each point (default 1 throughout); at a quadratic point the value is the norm over the two
conjugates.}
    if IsEmpty(Degrees) then Degrees := [1 : d in ds]; end if;
    require #Degrees eq #ds : "Degrees must have one entry per discriminant";
    Ldata := ShimuraCurveLattice(Xstar`D, Xstar`N);
    lambdas := ElementsOfNorm(Ldata`Q, [-d : d in ds], Ldata`O, Ldata`basis_L);
    table := [[* *] : f in fs];
    for j->d in ds do
        vals := AbsoluteValuesAtRationalCMPoint(fs, d, Xstar, Ldata : Lambda := lambdas[-d], PointDegree := Degrees[j]);
        for i->v in vals do
            Append(~table[i], RationalNumber(v));
        end for;
    end for;
    return table;
end intrinsic;

intrinsic CoordinateFromBorcherdsForms(vals0::List, vals1::List, ds::SeqEnum[RngIntElt],
                                       zeros0::SeqEnum[RngIntElt], zeros1::SeqEnum[RngIntElt]) -> List, FldRatElt, FldRatElt
{vals0, vals1 are the absolute values at the points ds (from BorcherdsFormValues) of two Borcherds
forms psi_0, psi_1 with the same double pole, psi_0 vanishing at zeros0 and psi_1 at zeros1.
Returns the values at ds of the function x = c_0 psi_0 with x = 0 at zeros0 and x = 1 at zeros1,
and the constants c_0, c_1 > 0 with |x| = c_0 |psi_0| and |1 - x| = c_1 |psi_1|.  The sign of x at
a point is the one for which |1 - x| agrees with c_1 |psi_1|; the call fails if neither or both do.}
    require #vals0 eq #ds and #vals1 eq #ds : "one value per discriminant";
    idx := func<d | Index(ds, d)>;
    require forall{d : d in zeros0 cat zeros1 | idx(d) ne 0} : "every zero must be among ds";
    require forall{d : d in zeros0 | vals0[idx(d)] eq 0} : "psi_0 does not vanish at zeros0";
    require forall{d : d in zeros1 | vals1[idx(d)] eq 0} : "psi_1 does not vanish at zeros1";
    // x = 1 at the zeros of psi_1 fixes c_0; 1 - x = 1 at the zeros of psi_0 fixes c_1
    c0s := {1/vals0[idx(d)] : d in zeros1};
    c1s := {1/vals1[idx(d)] : d in zeros0};
    require #c0s eq 1 : Sprintf("|psi_0| is not constant on the zeros of psi_1: 1/c_0 in %o", {1/c : c in c0s});
    require #c1s eq 1 : Sprintf("|psi_1| is not constant on the zeros of psi_0: 1/c_1 in %o", {1/c : c in c1s});
    c0 := Rep(c0s); c1 := Rep(c1s);
    xs := [* *];
    for j->d in ds do
        if is_infinite(vals0[j]) or is_infinite(vals1[j]) then
            require is_infinite(vals0[j]) and is_infinite(vals1[j]) :
                Sprintf("only one of psi_0, psi_1 has a pole at %o", d);
            Append(~xs, Infinity());
            continue;
        end if;
        a := c0*vals0[j]; b := c1*vals1[j];
        cands := {s*a : s in [1, -1] | AbsoluteValue(1 - s*a) eq b};
        require #cands eq 1 : Sprintf("sign of x at %o undetermined: |x| = %o, |1 - x| = %o", d, a, b);
        Append(~xs, Rep(cands));
    end for;
    return xs, c0, c1;
end intrinsic;

intrinsic BidegreeTwoRelation(xs::List, ys::List) -> RngMPolElt
{The relation F(x, y) = 0 of bidegree (2, 2) satisfied by the pairs (xs[i], ys[i]), where
Infinity() marks a pole, when the pairs determine it up to a scalar (nine coefficients, so at
least eight independent conditions).  F is returned with integer coefficients of content one and
positive leading coefficient.}
    require #xs eq #ys : "xs and ys must have the same length";
    rows := [];
    for i in [1..#xs] do
        X := is_infinite(xs[i]) select [Rationals() | 0, 1] else [Rationals() | 1, xs[i]];   // (X0 : X1), x = X1/X0
        Y := is_infinite(ys[i]) select [Rationals() | 0, 1] else [Rationals() | 1, ys[i]];
        Append(~rows, [X[1]^(2-a)*X[2]^a*Y[1]^(2-b)*Y[2]^b : a in [0..2], b in [0..2]]);
    end for;
    K := Kernel(Transpose(Matrix(Rationals(), rows)));
    require Dimension(K) eq 1 : Sprintf("the points leave a %o-dimensional space of relations", Dimension(K));
    kv := Basis(K)[1];
    P<x, y> := PolynomialRing(Rationals(), 2);
    mons := [x^a*y^b : a in [0..2], b in [0..2]];
    F := &+[kv[i]*mons[i] : i in [1..#mons]];
    F *:= LCM([Denominator(c) : c in Coefficients(F)]);
    F /:= GCD([Integers() | c : c in Coefficients(F)]);
    if LeadingCoefficient(F) lt 0 then F := -F; end if;
    return F;
end intrinsic;

intrinsic EllipticCurveFromRelation(F::RngMPolElt, xs::List, ys::List : BasePoint := 0) -> CrvEll, Map
{A minimal Weierstrass model of the genus-one plane curve F(x, y) = 0 and the map to it from the
projective closure of the curve, sending the base point to the origin.  The base point is the
BasePoint-th of the pairs (xs[i], ys[i]) or, by default, the first one that is finite and
nonsingular on the curve.}
    P := Parent(F);
    C := Curve(AffineSpace(P), F);
    require Genus(C) eq 1 : Sprintf("the relation defines a curve of genus %o", Genus(C));
    Cp := ProjectiveClosure(C);
    finite := [i : i in [1..#xs] | not is_infinite(xs[i]) and not is_infinite(ys[i])];
    if BasePoint ne 0 then
        require BasePoint in finite : "the base point must be finite";
        finite := [BasePoint];
    end if;
    require exists(pt){Cp![xs[i], ys[i], 1] : i in finite | IsNonsingular(Cp, Cp![xs[i], ys[i], 1])} :
        "no finite nonsingular point among the given pairs";
    E, phi := EllipticCurve(Cp, pt);
    Em, psi := MinimalModel(E);
    return Em, phi*psi;
end intrinsic;

intrinsic CMDivisorObstruction(Xstar::ShimuraQuot, ds::SeqEnum[RngIntElt] : PoleOrder := 0, Prec := 100) -> ModTupFld, RngIntElt
{For even D: the linear conditions a divisor sum_d a_d P_d, supported on the rational CM points of
discriminants ds of Xstar and written in point multiplicities, must satisfy to be the divisor of a
Borcherds form whose input has pole order at most PoleOrder at oo (default: the hauptmodul search's
window n_0 + k - 1 plus two powers of t).  Returns the space of such functionals, as row vectors
indexed like ds, and the pole order used.  The Borcherds forms here carry any weight, so the degree
is not among the conditions; on a genus-zero star the space is zero, and on a genus-one star whose
Jacobian has rank one it is spanned by one functional whose entry at a rational point is, up to a
common scale and shift, the position of the point in E(Q) modulo torsion (at X_0^142(1) the shift
is 2 and the positions are Guo-Yang's multiples of Q).  A divisor that fails a functional is not
principal; one that passes may still need a deeper pole order (the window grows with it).  The
discriminants may include quadratic points, whose entry is then the position of the whole cycle
with the shift counted twice, once per point.}
    require IsEven(Xstar`D) : "only the cusp-oo side is implemented (even D)";
    E, n, n0, t, eta_quotients := WeaklyHolomorphicBasis(Xstar`D, Xstar`N : Prec := Prec);
    k := -Valuation(qExpansionAtoo(t, 1));
    P := (PoleOrder eq 0) select n0 + k - 1 + 2*k else PoleOrder;
    ech_basis := basis_of_weakly_holomorphic_forms(P, eta_quotients, n0+1, n, t);
    mat, relevant_ds := coeffs_to_divisor_matrix(-P, Xstar`D, Xstar`N, Ncols(ech_basis));
    require forall{d : d in ds | Index(relevant_ds, -d) ne 0} :
        Sprintf("some of %o are out of reach at pole order %o (reachable: %o)", ds, P, [-d : d in relevant_ds]);
    ct := ChangeRing(ech_basis * ChangeRing(mat, BaseRing(ech_basis)), Rationals());
    V := RSpace(Rationals(), Ncols(ct));
    idx := [Index(relevant_ds, -d) : d in ds];
    W := sub<V | [V.j : j in idx]>;
    WI := W meet Image(ct);
    // the point multiplicity at an elliptic point of order e is the coefficient divided by e
    ells := NumberOfEllipticPointsByCMOrder(Xstar);
    e := [1 : d in ds];
    for q in Keys(ells) do
        for j->d in ds do
            if IsDefined(ells[q], d) then e[j] := q; end if;
        end for;
    end for;
    B := Matrix(Rationals(), #Basis(WI), #ds, [[(V!b)[idx[j]]/e[j] : j in [1..#ds]] : b in Basis(WI)]);
    return Kernel(Transpose(B)), P;
end intrinsic;
