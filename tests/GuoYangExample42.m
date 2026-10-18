// EXTERNAL SOURCE: Guo-Yang, "Equations of hyperelliptic Shimura curves", Compositio Math. 153
// (2017), Example 42, pp. 30-31 (journal version; the PDF is in Dropbox under
// Research/HyperellipticQuotients).  Everything compared against below is printed there: the
// principal parts of their forms f_1..f_4 in M^!(284), the four divisors, the table of x and y at
// ten CM points, the relation between x and y, the curve E142A1 and the positions of the CM
// points on it as multiples of Q = (1, 1).
//
// X = X_0^142(1)/W_{142,1} has genus one, so the hauptmodul search cannot succeed (no divisor of
// degree one is principal; this base sat in the "obstructed" list for that reason).  The paper
// builds the model from four Borcherds forms of degree two instead:
//     div psi_1 = P_-4  + P_-8  - 2 P_-3,     div psi_2 = P_-19 + P_-43 - 2 P_-3,
//     div psi_3 = P_-8  + P_-40 - 2 P_-20,    div psi_4 = P_-19 + P_-24 - 2 P_-20,
// x = 2^-10 psi_1 with 1 - x = psi_2, y = psi_3/2 with 1 - y = psi_4/2, and the relation
//     2 (x+1)^2 y^2 - (8x^2 + 11x + 1) y + 4x (2x+1) = 0,
// which is E142A1: y^2 + xy + y = x^3 - x^2 - 12x + 15.
// (The paper's sentence "y(tau_-8) = y(tau_-20) = 0" is a typo for y(tau_-8) = y(tau_-40) = 0; the
// divisor of psi_3 and the table both say so.)
//
// The forms found here differ from the paper's by forms with no CM divisor, which multiply the
// Borcherds form by a constant (measured: our psi_1, .., psi_4 are the paper's times 2^-18,
// 2^-10, 2^-3, 2^-1).  The values of x and y do not depend on this, and the constants are not
// checked against the paper, only printed and required to be powers of 2 like the paper's.

procedure test_guoyang_example42()
    printf "Guo-Yang Example 42: X_0^142(1)* (genus one) from Borcherds forms of degree two...";
    D := 142; N := 1;
    Xstar := CreateShimuraQuot(D, N, Set(Divisors(D*N)));
    Xstar`g := GenusShimuraCurveQuotient(D, N, Xstar`W); Xstar`CurveID := 0;
    assert Xstar`g eq 1;
    curves := GetQuotientsAndGenera([Xstar]);
    _ := exists(star){c : c in curves | IsStarCurve(c)};

    // the ten rational CM points of the paper's table, and no others up to class number 2
    ds := [-3, -4, -8, -19, -20, -24, -40, -43, -148, -232];
    rat, _ := RationalandQuadraticCMPoints(star : bd := 2, coprime_to_level := false);
    assert Sort([p[1] : p in rat]) eq Sort(ds);

    paper_div := [{<-4, 1>, <-8, 1>, <-3, -2>}, {<-19, 1>, <-43, 1>, <-3, -2>},
                  {<-8, 1>, <-40, 1>, <-20, -2>}, {<-19, 1>, <-24, 1>, <-20, -2>}];

    // (1) The paper's own forms, as printed (principal parts at oo; the constant terms are 0),
    // have the printed divisors under this repo's reading of a principal part, which divides the
    // coefficient at an elliptic point by its order and collects every exponent m r^2 that meets
    // the same discriminant.  This checks the convention before any form of ours is trusted.
    R<q> := LaurentSeriesRing(Rationals());
    paper_pp := [
        [<-87,-2>, <-71,-2>, <-48,-2>, <-36,-2>, <-16,2>, <-15,-2>, <-12,-2>, <-9,2>, <-7,-2>, <-3,-2>, <-2,2>, <-1,2>],
        [<-116,2>, <-87,-1>, <-79,-1>, <-71,2>, <-60,2>, <-48,-2>, <-43,1>, <-29,-2>, <-19,1>, <-15,-4>, <-12,-2>, <-7,2>, <-3,-2>],
        [<-87,1>, <-79,-2>, <-76,1>, <-71,-2>, <-48,2>, <-40,-1>, <-32,3>, <-20,-2>, <-19,-1>, <-12,-2>, <-10,2>, <-8,1>, <-7,-2>, <-2,-2>],
        [<-79,-1>, <-76,1>, <-48,1>, <-40,-1>, <-32,1>, <-20,-1>, <-12,-1>, <-10,1>, <-6,1>, <-5,-1>, <-2,-1>]];
    for i in [1..4] do
        f := &+[R | c[2]*q^c[1] : c in paper_pp[i]] + O(q);
        assert Set(DivisorOfBorcherdsForm(f, star)) eq paper_div[i];
    end for;

    // (2) Forms with the paper's divisors from the prescribed-divisor search.  Their pole orders
    // at oo (87, 116, 87, 79 in the paper) exceed the hauptmodul search's window of 54, which is
    // what the pole-order ladder is for.  The principal parts must be integral, so that each
    // Borcherds form is a function and not a root of one.
    divisors := [<[-4, -8], -3>, <[-19, -43], -3>, <[-8, -40], -20>, <[-19, -24], -20>];
    fs := BorcherdsForms(star, curves : Divisors := divisors);
    assert Keys(fs) eq {1, 2, 3, 4};
    for i in [1..4] do
        assert Set(DivisorOfBorcherdsForm(fs[i], star)) eq paper_div[i];
        g := qExpansionAtoo(fs[i], 1);
        assert forall{m : m in [Valuation(g)..0] | IsIntegral(Coefficient(g, m))};
        assert Coefficient(g, 0) eq 0;
    end for;

    // (3) Schofer values at the ten points, the two coordinates, and the paper's table entry by
    // entry.  Every disagreement is printed before the assertion fires.
    vals := BorcherdsFormValues([fs[i] : i in [1..4]], ds, star);
    xs, cx0, cx1 := CoordinateFromBorcherdsForms(vals[1], vals[2], ds, [-4, -8], [-19, -43]);
    ys, cy0, cy1 := CoordinateFromBorcherdsForms(vals[3], vals[4], ds, [-8, -40], [-19, -24]);
    paper_x := [* Infinity(), 0, 0, 1, -1, 1/2, -1/2, 1, -1, -1/2 *];
    paper_y := [* 2, 1/2, 0, 1, Infinity(), 1, 0, 3/2, -2, -5 *];
    same := func<a, b | (Type(a) eq Infty and Type(b) eq Infty) or
                        (Type(a) ne Infty and Type(b) ne Infty and a eq b)>;
    bad := [* *];
    for j->d in ds do
        if not same(xs[j], paper_x[j]) then Append(~bad, <"x", d, xs[j], paper_x[j]>); end if;
        if not same(ys[j], paper_y[j]) then Append(~bad, <"y", d, ys[j], paper_y[j]>); end if;
    end for;
    if #bad gt 0 then
        printf "\n  entries disagreeing with the paper's table (coordinate, d, ours, paper): %o\n", bad;
    end if;
    assert #bad eq 0;
    // the normalisation constants of our forms (the paper's are 2^-10, 1, 1/2, 1/2; see the header)
    printf " [x = %o psi_1, 1 - x = %o psi_2, y = %o psi_3, 1 - y = %o psi_4]", cx0, cx1, cy0, cy1;
    assert forall{c : c in [cx0, cx1, cy0, cy1] | Set(PrimeDivisors(Numerator(c)*Denominator(c))) subset {2}};

    // (4) The relation between x and y, and the curve.
    F := BidegreeTwoRelation(xs, ys);
    P<x, y> := Parent(F);
    assert F eq 2*(x+1)^2*y^2 - (8*x^2 + 11*x + 1)*y + 4*x*(2*x+1);
    // base point P_-8 = (0, 0), which the paper's table sends to the origin
    E, toE := EllipticCurveFromRelation(F, xs, ys : BasePoint := Index(ds, -8));
    assert CremonaReference(E) eq "142a1";
    E142A1 := EllipticCurve([1, -1, 1, -12, 15]);
    assert IsIsomorphic(E, E142A1);

    // (5) The CM points on E142A1: the paper lists them as -Q, -2Q, O, Q, 2Q, 3Q, 4Q, -3Q, -4Q,
    // -6Q with Q = (1, 1).  Both maps send P_-8 to O, so they differ at most by [-1].  The eight
    // finite points go through the map; the two poles P_-3 and P_-20 sit at infinity of the plane
    // model, where it is singular, so they are placed by the divisors of x and y instead:
    // P_-4 + P_-8 - 2 P_-3 and P_-8 + P_-40 - 2 P_-20 are principal, and E(Q) is torsion-free
    // (142a1 has trivial torsion), so each pole is the unique half of the sum of the two zeros.
    Q := E![1, 1, 1];
    assert #TorsionSubgroup(E) eq 1;
    paper_mult := [-1, -2, 0, 1, 2, 3, 4, -3, -4, -6];
    multiple := function(Pt)
        for m in [-8..8] do if m*Q eq Pt then return m; end if; end for;
        return 10^6;    // not a small multiple of Q
    end function;
    Cp := Domain(toE);
    mults := AssociativeArray();
    for j->d in ds do
        if Type(xs[j]) eq Infty or Type(ys[j]) eq Infty then continue; end if;
        mults[d] := multiple(toE(Cp![xs[j], ys[j], 1]));
    end for;
    assert IsEven(mults[-4] + mults[-8]) and IsEven(mults[-8] + mults[-40]);
    mults[-3] := (mults[-4] + mults[-8]) div 2;
    mults[-20] := (mults[-8] + mults[-40]) div 2;
    ours := [mults[d] : d in ds];
    assert ours eq paper_mult or ours eq [-m : m in paper_mult];

    // (6) Without any Schofer value: the one linear condition a CM divisor must satisfy at pole
    // order 102 (the window in which the paper's forms live) is the Abel-Jacobi condition, and
    // its coefficients are the positions above up to a common shift and scale.
    K, P := CMDivisorObstruction(star, ds : PoleOrder := 102);
    assert Dimension(K) eq 1;
    v := Eltseq(Basis(K)[1]);
    ones := [1 : d in ds];
    assert Rank(Matrix(Rationals(), [v, paper_mult, ones])) eq 2;
    assert Rank(Matrix(Rationals(), [v, ones])) eq 2;
    // and each of the paper's divisors passes it, read in point multiplicities
    for dv in paper_div do
        assert &+[c[2]*v[Index(ds, c[1])] : c in dv] eq 0;
    end for;
    printf "Done!\n";
end procedure;

test_guoyang_example42();
