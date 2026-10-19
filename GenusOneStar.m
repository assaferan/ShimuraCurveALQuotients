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
// Where Guo-Yang's Assumption 40 enters, and when it is automatic.  A Borcherds form transforms under the Atkin-Lehner
// group by a character, and it is a function on the star quotient only when that character is
// trivial, which Guo-Yang assume (their Assumption 40; Example 42 is stated under it).  The
// search here does not check the character: it only finds a form whose divisor on the star is the
// prescribed one.  Some power of the form is a function, so the prescribed divisor, which is
// rational, has a class of finite order in E(Q).  When E(Q) has no torsion the class is trivial,
// a function g with that divisor exists, and the form divided by g has no zeros or poles, so the
// character is trivial and the assumption holds; this is the case at X_0^142(1)* (142a1 has
// trivial torsion).  When E(Q) has torsion the form may be a root of a function, its Schofer
// values are then those of the root, and a model built from them is conditional on the
// assumption.  The character can be dealt with exactly.  W is an elementary abelian 2-group, so
// every character has order at most 2 and the square of every form is a function.  A form with
// character chi != 1 is a function on X/ker(chi), a quotient by an index-2 subgroup of W, and its
// square has even divisor on the star, so X/ker(chi) is an etale double cover of the star and
// has genus one.  Hence: when no index-2 quotient of the star has genus one, every Borcherds form
// is a function, whatever the torsion of E(Q) (EtaleDoubleCovers).  When such a quotient X' exists,
// the form psi with divisor A - P (A - P of order 2 in E(Q)) is an anti-invariant function on X',
// its square is the function with divisor 2A - 2P in L(2P), and the sign of that square at a
// rational CM point Q of the star is the sign of the discriminant of the field of definition of
// the preimages of Q on X' (a rational preimage gives a rational value of psi, a quadratic pair
// gives r*sqrt(delta) and its square r^2 delta), up to one global sign (SquareCoordinateSigns).
// The coordinates built this way use only functions, so the model does not depend on the
// assumption.

import "BorcherdsForms.m" : basis_of_weakly_holomorphic_forms, coeffs_to_divisor_matrix;

function is_infinite(v)
    return Type(v) eq Infty;
end function;

intrinsic BorcherdsFormValues(fs::SeqEnum[EtaQuot], ds::SeqEnum[RngIntElt], Xstar::ShimuraQuot : Degrees := []) -> SeqEnum
{The absolute values |psi_f| of the Borcherds forms of the forms in fs at the CM points of Xstar of
discriminants ds, one list per form (a zero gives 0, a pole Infinity()); a value that is not a
rational number is an error.  Degrees gives the degree of each point (default 1 throughout); at
a quadratic point the value is the norm over the two conjugates.}
    if IsEmpty(Degrees) then Degrees := [1 : d in ds]; end if;
    require #Degrees eq #ds : "Degrees must have one entry per discriminant";
    Ldata := ShimuraCurveLattice(Xstar`D, Xstar`N);
    lambdas := ElementsOfNorm(Ldata`Q, [-d : d in ds], Ldata`O, Ldata`basis_L);
    table := [[* *] : f in fs];
    for j->d in ds do
        vals := AbsoluteValuesAtRationalCMPoint(fs, d, Xstar, Ldata : Lambda := lambdas[-d], PointDegree := Degrees[j]);
        for i->v in vals do
            // a fractional exponent means the form is not a function on the star (a nontrivial
            // character) or the local factor of Schofer's formula at a prime dividing both d and
            // N is wrong; either way nothing downstream is valid, so stop here
            ok := true;
            try r := RationalNumber(v); catch e ok := false; end try;
            error if not ok, Sprintf("BorcherdsFormValues: the value of form %o at discriminant %o is not a rational number: %o", i, d, v);
            Append(~table[i], r);
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

// ---- certification of a genus-three double cover of the star ----------------------------
intrinsic BiellipticBranchCriterion(B::SeqEnum, L::PtEll) -> BoolElt, List
{B is the branch locus of a double cover C -> E of a genus-one curve by a genus-three curve: four
distinct points of E over one number field.  L is a point of E giving the square-root class of
the cover: the cover is E(sqrt phi) with div phi = B - 2 D and L = the sum of the two points of
D (for a Borcherds form with divisor B - 4 P_oo, L = 2 P_oo).  Returns whether C is
hyperelliptic, and the list of pairings of B into two pairs with equal sums, each with its sum Q0.
C is hyperelliptic exactly when some pairing has Q0 = L.  If C is hyperelliptic its hyperelliptic
involution commutes with the covering involution and descends to tau: P -> Q0 - P on E, which
preserves B; tau fixes no branch point, since over a fixed branch point a lift of tau to C would
fix at most one point and have at most seven fixed points in all instead of eight; so tau swaps
B in two pairs, each of sum Q0.  Conversely the four double covers branched at B have classes
Q0 + E[2], and only the class Q0, the pull-back of O(1) from E/tau = P^1, carries a lift of tau
with eight fixed points; the other three are Ciani quartics.  When no pairing has equal sums the
verdict needs no class at all.}
    require #B eq 4 and #Set(B) eq 4 : "four distinct branch points";
    PS := Universe(B);
    LK := PS!Eltseq(L);
    pairings := [[[1,2],[3,4]], [[1,3],[2,4]], [[1,4],[2,3]]];
    equal := [* *];
    for pr in pairings do
        s1 := B[pr[1][1]] + B[pr[1][2]];
        s2 := B[pr[2][1]] + B[pr[2][2]];
        if s1 eq s2 then Append(~equal, <pr, s1>); end if;
    end for;
    return exists{e : e in equal | e[2] eq LK}, equal;
end intrinsic;

intrinsic IsHyperellipticByCanonicalQuadric(FC::FldFun) -> BoolElt
{Whether the genus-three function field FC is hyperelliptic: its three holomorphic differentials
satisfy a quadratic relation (the canonical image is a conic) exactly then.  Decided by a
constant-coefficient linear relation among the six pairwise products, written as functions by
dividing by one differential; the products are expanded in a basis of FC over the rational
function field, which turns the relation into linear algebra over the constants.}
    require Genus(FC) eq 3 : "genus three only";
    V, m := SpaceOfHolomorphicDifferentials(FC);
    omegas := [m(V.i) : i in [1..3]];
    fs := [omegas[i]/omegas[1] : i in [1..3]];       // functions: omega_i = f_i omega_1
    prods := [fs[i]*fs[j] : j in [i..3], i in [1..3]];
    // coordinates over the rational function field at the bottom of the tower
    flatten := function(f)
        cs := [f];
        while not ISA(Type(Universe(cs)), FldFunRat) do
            cs := &cat[Eltseq(c) : c in cs];
        end while;
        return cs;
    end function;
    vecs := [flatten(f) : f in prods];
    den := LCM(&cat[[Denominator(c) : c in v] : v in vecs]);
    polys := [[Numerator(den*c) : c in v] : v in vecs];
    degs := [Maximum([Degree(q) : q in polys[i]]) : i in [1..#polys]];
    dmax := Maximum(degs cat [0]);
    rows := [&cat[[Coefficient(q, e) : e in [0..dmax]] : q in polys[i]] : i in [1..#polys]];
    kern := Kernel(Matrix(BaseRing(Universe(vecs[1])), rows));
    for v in Basis(kern) do
        if IsZero(&+[v[i]*prods[i] : i in [1..6]]) then return true; end if;
    end for;
    return false;
end intrinsic;

// ---- characters: the etale double covers of the star, and signs of square-coordinates ----------

intrinsic EtaleDoubleCovers(Xstar::ShimuraQuot, curves::SeqEnum[ShimuraQuot]) -> SeqEnum
{The quotients of X_0(D, N) by index-2 subgroups of the Atkin-Lehner group of the genus-one star
Xstar that have genus one, i.e. the etale double covers of the star among the Atkin-Lehner
quotients.  A Borcherds form with nontrivial character is a function on one of these; when the
list is empty every Borcherds form on Xstar is a function.}
    require Xstar`g eq 1 : "the star must have genus one";
    return [curves[i] : i in Xstar`CoveredBy | curves[i]`g eq 1];
end intrinsic;

intrinsic SquareCoordinateSigns(Xp::ShimuraQuot, ds::SeqEnum[RngIntElt]) -> SeqEnum
{For the function x = psi^2 on the star, psi an anti-invariant function on the etale double cover
Xp of the star: the sign of x at the rational CM point of each discriminant in ds, up to one
common sign.  It is +1 when the preimages of the point on Xp are rational, and the sign of the
discriminant of their quadratic field of definition otherwise.}
    signs := [];
    for d in ds do
        deg := DegreeOfFieldOfDefinitionOfCMPoint(Xp, d);
        require deg in {1, 2} : Sprintf("the CM points of discriminant %o on the cover have degree %o", d, deg);
        if deg eq 1 then
            Append(~signs, 1);
        else
            K := FieldOfDefinitionOfCMPoint(Xp, d);
            require Degree(K) eq 2 : Sprintf("the field of definition at %o is not quadratic", d);
            Append(~signs, Sign(Discriminant(MaximalOrder(K))));
        end if;
    end for;
    return signs;
end intrinsic;
