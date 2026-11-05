// HyperellipticDetection.m
//
// Whether a curve given by its function field is hyperelliptic, and a model y^2 = f(x) for a
// fibre product of double covers of the line that is a double cover of a conic with a rational
// point.
//
// The test is Max Noether's: the canonical map sends a curve of genus g >= 3 to P^(g-1); a
// hyperelliptic curve goes two-to-one onto a rational normal curve of degree g-1, which lies on
// (g-1)(g-2)/2 independent quadrics, while the canonical image of a non-hyperelliptic curve lies
// on exactly (g-2)(g-3)/2.  So the number of linear relations among the pairwise products of a
// basis of holomorphic differentials decides, and nothing else is needed: no model of the curve
// beyond its function field, no search for an involution.

intrinsic IsHyperellipticFunctionField(K::FldFun) -> BoolElt, RngIntElt
{Whether the curve with function field K (genus at least 2) is hyperelliptic, by the number of
 quadratic relations among its holomorphic differentials: (g-1)(g-2)/2 if it is, (g-2)(g-3)/2 if
 not.  Also returns that number.}
    g := Genus(K);
    require g ge 2 : "the curve must have genus at least 2";
    if g eq 2 then return true, 0; end if;
    V, mp := SpaceOfHolomorphicDifferentials(K);
    ws := [mp(V.i) : i in [1..g]];
    phis := [ws[i] / ws[1] : i in [1..g]];
    prods := [phis[i]*phis[j] : j in [i..g], i in [1..g]];
    nrel := Dimension(Relations(prods, BaseRing(V)));
    hyp := nrel eq (g-1)*(g-2) div 2;
    require hyp or nrel eq (g-2)*(g-3) div 2 : "the number of quadratic relations matches neither case";
    return hyp, nrel;
end intrinsic;

// The squarefree part of a polynomial together with its leading constant: f = c * g * h^2 with
// g monic squarefree; returns c*g.  (SquarefreePart drops the constant, which decides the twist.)
function squarefree_with_constant(f)
    fac, u := Factorization(f);
    return u * &*[Parent(f) | fa[1]^(fa[2] mod 2) : fa in fac];
end function;

intrinsic ConicFactorOfFibreProduct(fs::SeqEnum[RngUPolElt]) -> BoolElt, RngUPolElt, SetEnum
{For the fibre product over the t-line of the double covers y_i^2 = f_i(t): a nonempty subset S
 of the factors whose product is, up to squares, a polynomial c of degree at most 2, so that
 y^2 = c(t) is a genus-0 quotient of the curve by one of its (Z/2)^k-automorphisms.  Returns
 false when there is none; otherwise true, c and S.}
    for S in Subsets({1..#fs}) do
        if IsEmpty(S) then continue; end if;
        c := squarefree_with_constant(&*[fs[i] : i in S]);
        if Degree(c) le 2 and Degree(c) ge 1 then return true, c, S; end if;
    end for;
    return false, _, _;
end intrinsic;

intrinsic ConicHasRationalPoint(c::RngUPolElt) -> BoolElt
{Whether the genus-0 curve y^2 = c(t), deg c <= 2, has a rational point.}
    require Degree(c) le 2 : "c must have degree at most 2";
    if Degree(c) le 1 then return true; end if;
    P2<T, Y, Z> := ProjectiveSpace(Rationals(), 2);
    return HasRationalPoint(Conic(P2, Y^2 - &+[Coefficient(c, i) * T^i * Z^(2 - i) : i in [0..2]]));
end intrinsic;

// Small integral models: x -> (x + b)/a, y -> y a^(d/2) over a few a and b, the square part of the
// content divided out; the one with the shortest coefficients is kept.
function small_model(f)
    P := Parent(f); x := P.1;
    d := Degree(f); dd := IsEven(d) select d else d + 1;
    best := f; bestsize := &+[#Sprint(c) : c in Coefficients(f)];
    for a in [1, 2, 3, 4, 5, 6, 8, 10, 12, 20, 25, 50, 125] do
        for b in [0, 1, -1, 2, -2] do
            g := Evaluate(f, (x + b)/a) * a^dd;
            L := LCM([Denominator(c) : c in Coefficients(g)]);
            g := g * L^2;
            cont := GCD([Integers() | c : c in Coefficients(g)]);
            sq := &*[Integers() | fa[1]^(fa[2] div 2) : fa in Factorization(cont)];
            g := g / sq^2;
            size := &+[#Sprint(c) : c in Coefficients(g)];
            if size lt bestsize then best := g; bestsize := size; end if;
        end for;
    end for;
    return best;
end function;

intrinsic HyperellipticModelOverConic(fs::SeqEnum[RngUPolElt]) -> CrvHyp
{For the fibre product over the t-line of the double covers y_i^2 = f_i(t), when it is a double
 cover of a genus-0 quotient y^2 = c(t) with a rational point (ConicFactorOfFibreProduct and
 ConicHasRationalPoint): a model y^2 = f(x) over Q.  The conic is parametrised by P^1 from a
 rational point, t = t(x), and the remaining factor y_1^2 = f_1(t(x)) is cleared of square
 denominators; the coefficients are then made integral and small by a change of x.  Only the
 case of two factors is implemented, which is the case of a curve of the form y_1^2 = f_1(t),
 y_2^2 = c(t).}
    require #fs eq 2 : "two factors only";
    ok, c, S := ConicFactorOfFibreProduct(fs);
    require ok : "no factor is a conic";
    require ConicHasRationalPoint(c) : "the conic has no rational point";
    i_c := Representative(S);
    require #S eq 1 : "the conic must be one of the two factors";
    f1 := fs[3 - i_c];
    F<X> := FunctionField(Rationals());
    if Degree(c) le 1 then
        // y^2 = c(t) is a line: t is a rational function of y
        tx := Degree(c) eq 0 select X else (X^2 - Coefficient(c, 0)) / Coefficient(c, 1);
    else
        P2<T, Y, Z> := ProjectiveSpace(Rationals(), 2);
        C0 := Conic(P2, Y^2 - &+[Coefficient(c, i) * T^i * Z^(2 - i) : i in [0..2]]);
        _, pt := HasRationalPoint(C0);
        par := Parametrization(C0, pt);
        im := [Evaluate(d, [X, 1]) : d in DefiningPolynomials(par)];
        tx := im[1] / im[3];
    end if;
    fx := Evaluate(f1, tx);
    g := Numerator(fx) * Denominator(fx);                    // y'^2 = g with y' = y_1 * denominator
    f := squarefree_with_constant(PolynomialRing(Rationals())!g);
    return HyperellipticCurve(small_model(f));
end intrinsic;
