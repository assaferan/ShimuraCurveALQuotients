// X_0^21(58)* IS HYPERELLIPTIC, PART ONE: THE POINTS, THE JACOBIAN, THE TORSION.
//
// The star X* = X_0^21(58)* has genus three.  Its rational CM points of class number at most four
// include the eight points of discriminants -16, -120, -88, -100, -112, -372, -168, -228, which a
// numerical screen pairs into four fibres F2, F4, F5, F6 of the putative hyperelliptic map.  The
// offline test tests/_offline/Hyperelliptic21_58.m certifies that the differences of these fibres are
// torsion in Jac(X*)(Q), through the lemma of paper/hyperelliptic-21-58.tex, which needs:
//   (a) each of the eight discriminants carries exactly one point of X*, and it is rational
//       (the CM points of one discriminant form one Galois orbit, so this is the degree 1 of their
//       field of definition);
//   (b) Jac(X*) is isogenous to A_f for a newform f of level 609 = 21 * 29, i.e. old at 2 (the
//       correspondence to X_0^21(29)* is then injective on Jac(X*) (x) Q); f is identified among the
//       finitely many newform orbits of dimension three and level dividing 1218 that Jacquet-Langlands
//       allows by its Hecke polynomials at 5, 11, 13;
//   (c) the torsion of Jac(X*)(Q) has order at most 2, so four effective divisors of degree two in
//       one torsion coset contain two that are linearly equivalent.
// This test checks (a), (b) and (c), from the trace formula and modular symbols, with no model.

procedure test_hyperelliptic_21_58()
    printf "X_0^21(58)*: eight rational CM points, Jacobian A_f with f = 609.2.a.d, torsion bound 2...";
    D := 21; N := 58;
    W := {m : m in Divisors(D*N) | GCD(m, (D*N) div m) eq 1};
    Xstar := CreateShimuraQuot(D, N, W);
    Xstar`g := GenusShimuraCurveQuotient(D, N, W); Xstar`CurveID := 0;
    assert Xstar`g eq 3;
    curves := GetQuotientsAndGenera([Xstar]);
    _ := exists(star){c : c in curves | IsStarCurve(c)};

    // (a) the eight points
    // (the class-number table behind RationalandQuadraticCMPoints lists fundamental discriminants
    // only; -16, -100, -112 are not fundamental, so the degree is computed directly)
    pts := [-16, -120, -88, -100, -112, -372, -168, -228];
    assert forall{d : d in pts | DegreeOfFieldOfDefinitionOfCMPoint(star, d) eq 1};

    // (b) the Hecke polynomials of Jac(X*) at 5, 11, 13, from the point counts over F_p, F_{p^2}, F_{p^3}
    R<x> := PolynomialRing(Rationals());
    g := 3;
    function HeckePoly(p)
        s := [ p^i + 1 - ComputePointsViaTrace(Xstar, p, i) : i in [1..g] ];
        e := [Rationals() | 1];
        for k in [1..g] do Append(~e, (&+[ (-1)^(i-1) * e[k-i+1] * s[i] : i in [1..k] ]) / k); end for;
        E := func< k | k le g select e[k+1] else p^(k-g) * e[2*g-k+1] >;
        P := &+[ (-1)^k * E(k) * x^(2*g-k) : k in [0..2*g] ];
        // P(x) = x^g Q(x + p/x): solve for Q
        Mat := ZeroMatrix(Rationals(), g+1, 2*g+1);
        for k in [0..g] do for j in [0..k] do Mat[k+1, g + k - 2*j + 1] +:= Binomial(k, j) * p^j; end for; end for;
        ok, q := IsConsistent(Mat, Vector(Rationals(), [ Coefficient(P, i) : i in [0..2*g] ]));
        assert ok;
        return R![ q[k+1] : k in [0..g] ];
    end function;
    polys := [ HeckePoly(p) : p in [5, 11, 13] ];
    assert polys eq [ x^3 + 6*x^2 + 8*x + 2, x^3 + 10*x^2 + 28*x + 20, x^3 + 6*x^2 + 12*x + 8 ];
    assert IsIrreducible(polys[1]);   // so Jac(X*) is isotypic for one orbit of dimension three
    // the candidates: newform orbits of weight 2 and dimension 3 at the levels 21 * m, m | 58, with the
    // Atkin-Lehner signs that make them W-invariant on the Shimura curve (the sign flips at 3 and 7)
    matches := [];
    for m in Divisors(58) do
        lev := 21 * m;
        Snew := NewSubspace(CuspidalSubspace(ModularSymbols(lev, 2, +1)));
        for A in NewformDecomposition(Snew) do
            if Dimension(A) ne 3 then continue; end if;
            if [ HeckePolynomial(A, p) : p in [5, 11, 13] ] ne polys then continue; end if;
            signs := [ <q, Integers()!AtkinLehnerOperator(A, q)[1, 1]> : q in PrimeDivisors(lev) ];
            Append(~matches, <lev, signs>);
        end for;
    end for;
    assert #matches eq 1 and matches[1][1] eq 609;
    assert Set(matches[1][2]) eq { <3, -1>, <7, -1>, <29, 1> };

    // (c) the torsion bound from the trace formula
    t := 0;
    for p in [q : q in PrimesInInterval(5, 60) | (D*N) mod q ne 0] do
        s := [ p^i + 1 - ComputePointsViaTrace(Xstar, p, i) : i in [1..g] ];
        e := [Rationals() | 1];
        for k in [1..g] do Append(~e, (&+[ (-1)^(i-1) * e[k-i+1] * s[i] : i in [1..k] ]) / k); end for;
        E := func< k | k le g select e[k+1] else p^(k-g) * e[2*g-k+1] >;
        t := GCD(t, Integers()!(&+[ (-1)^k * E(k) : k in [0..2*g] ]));
    end for;
    assert t eq 2;
    printf " [torsion bound %o, four fibres] Done!\n", t;
end procedure;

test_hyperelliptic_21_58();
