// X_0^39(34)* IS HYPERELLIPTIC, PART ONE: THE POINTS, THE JACOBIAN, THE TORSION.
//
// The star X* = X_0^39(34)* has genus three.  Five effective divisors of degree two supported on CM
// points -- the pairs of rational points {P_-84, P_-132}, {P_-15, P_-52}, {P_-1092, P_-372} and the
// quadratic points Z_-3432, Z_-1560 (the two conjugate CM points of those discriminants) -- have
// pairwise torsion differences in Jac(X*)(Q), certified by tests/_offline/Hyperelliptic39_34.m
// through the lemma of paper/hyperelliptic-21-58.tex (for the odd discriminant -15 through the
// plus-type dictionary of the level-17 curve, paper/hyperelliptic-39-34.tex).  The lemma and the
// counting need:
//   (a) the six rational discriminants carry one rational point each, and -3432, -1560 one Galois
//       orbit of two points each (the degree of the field of definition);
//   (b) Jac(X*) is isogenous to A_f for a newform f of level 663 = 39 * 17, i.e. old at 2, identified
//       among the newform orbits of dimension three at the levels 39m, m | 34, by its Hecke
//       polynomials at 7, 11, 19;
//   (c) the torsion of Jac(X*)(Q) has order at most 4, so five effective divisors of degree two in one
//       torsion coset contain two that are linearly equivalent, and the curve has a g^1_2.

procedure test_hyperelliptic_39_34()
    printf "X_0^39(34)*: the CM points, Jacobian A_f with f of level 663, torsion bound 4...";
    D := 39; N := 34;
    W := {m : m in Divisors(D*N) | GCD(m, (D*N) div m) eq 1};
    Xstar := CreateShimuraQuot(D, N, W);
    Xstar`g := GenusShimuraCurveQuotient(D, N, W); Xstar`CurveID := 0;
    assert Xstar`g eq 3;
    curves := GetQuotientsAndGenera([Xstar]);
    _ := exists(star){c : c in curves | IsStarCurve(c)};

    // (a) the points
    assert forall{d : d in [-84, -132, -15, -52, -1092, -372] | DegreeOfFieldOfDefinitionOfCMPoint(star, d) eq 1};
    assert forall{d : d in [-3432, -1560] | DegreeOfFieldOfDefinitionOfCMPoint(star, d) eq 2};

    // (b) the Hecke polynomials at 7, 11, 19 from the point counts over F_p, F_{p^2}, F_{p^3}
    R<x> := PolynomialRing(Rationals());
    g := 3;
    function Frob(p)
        s := [ p^i + 1 - ComputePointsViaTrace(Xstar, p, i) : i in [1..g] ];
        e := [Rationals() | 1];
        for k in [1..g] do Append(~e, (&+[ (-1)^(i-1) * e[k-i+1] * s[i] : i in [1..k] ]) / k); end for;
        return func< k | k le g select e[k+1] else p^(k-g) * e[2*g-k+1] >;
    end function;
    function HeckePoly(p)
        E := Frob(p);
        P := &+[ (-1)^k * E(k) * x^(2*g-k) : k in [0..2*g] ];
        Mat := ZeroMatrix(Rationals(), g+1, 2*g+1);
        for k in [0..g] do for j in [0..k] do Mat[k+1, g + k - 2*j + 1] +:= Binomial(k, j) * p^j; end for; end for;
        ok, q := IsConsistent(Mat, Vector(Rationals(), [ Coefficient(P, i) : i in [0..2*g] ]));
        assert ok;
        return R![ q[k+1] : k in [0..g] ];
    end function;
    polys := [ HeckePoly(p) : p in [7, 11, 19] ];
    assert polys eq [ x^3 + 2*x^2 - 4*x - 4, x^3 + 8*x^2 + 8*x - 16, x^3 + 4*x^2 - 8*x - 16 ];
    assert IsIrreducible(polys[1]);   // so Jac(X*) is isotypic for one orbit of dimension three
    matches := [];
    for m in Divisors(34) do
        lev := 39 * m;
        Snew := NewSubspace(CuspidalSubspace(ModularSymbols(lev, 2, +1)));
        for A in NewformDecomposition(Snew) do
            if Dimension(A) ne 3 then continue; end if;
            if [ HeckePolynomial(A, p) : p in [7, 11, 19] ] ne polys then continue; end if;
            Append(~matches, <lev, [ <q, Integers()!AtkinLehnerOperator(A, q)[1, 1]> : q in PrimeDivisors(lev) ]>);
        end for;
    end for;
    assert #matches eq 1 and matches[1][1] eq 663;
    assert Set(matches[1][2]) eq { <3, -1>, <13, -1>, <17, 1> };

    // (c) the torsion bound from the trace formula
    t := 0;
    for p in [q : q in PrimesInInterval(5, 60) | (D*N) mod q ne 0] do
        E := Frob(p);
        t := GCD(t, Integers()!(&+[ (-1)^k * E(k) : k in [0..2*g] ]));
    end for;
    assert t eq 4;
    printf " [torsion bound %o, five divisors] Done!\n", t;
end procedure;

test_hyperelliptic_39_34();
