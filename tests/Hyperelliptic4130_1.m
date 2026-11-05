// X_0^4130(1)* IS HYPERELLIPTIC, PART ONE: THE STAR, ITS RATIONAL CM POINTS, AND THE TORSION OF ITS JACOBIAN.
// About four minutes.
//
// Y = X_0^4130(1)* = X_0^4130(1)/W (W the full Atkin-Lehner group, order 16) has genus three and no further quotient,
// so it is hyperelliptic exactly when it carries two distinct linearly equivalent divisors of degree two.  Four
// candidates come from the Shimizu screen of the invariant differentials at the rational CM points (not used by the
// proof): D_1 = P_-35 + P_-532, D_2 = P_-163 + P_-595, D_3 = P_-315 + P_-520, D_4 = P_-340 + P_-420.
// tests/_offline/Hyperelliptic4130_1.m certifies, by the dual test of paper/hyperelliptic-15-38.tex, Section 2, that
// their differences are TORSION in J(Y)(Q).  This test shows that is enough (paper/hyperelliptic-4130-1.tex): the
// rational torsion of J(Y) is trivial, since the orders #J(Y)(F_p) = L_p(1), from the point counts over F_p, F_{p^2},
// F_{p^3} by the trace formula, have gcd 1 over the good primes below 120.  So a torsion difference is principal, and
// two distinct degree-two divisors in one class give a g^1_2.  The test also checks that the eight CM points are
// rational on Y and records the traces of Frobenius at 3 used by the offline test to recognise the space W: T(9)
// has characteristic polynomial x^3 + 3x^2 - 3 on it.

procedure test_hyperelliptic_4130_1()
    printf "X_0^4130(1)*: the eight CM points are rational and J(Y)(Q) is torsion-free, so one certified fibre difference is a g^1_2...";
    D := 4130; N := 1;
    W := {m : m in Divisors(D*N) | GCD(m, (D*N) div m) eq 1};
    assert #W eq 16;
    Y := CreateShimuraQuot(D, N, W); Y`g := GenusShimuraCurveQuotient(D, N, W); Y`CurveID := 0;
    assert Y`g eq 3;
    curves := GetQuotientsAndGenera([Y]);
    _ := exists(star){c : c in curves | IsStarCurve(c)};
    discs := [-35, -532, -163, -595, -315, -520, -340, -420];
    assert forall{d : d in discs | DegreeOfFieldOfDefinitionOfCMPoint(star, d) eq 1};
    assert -315 eq 9 * (-35) and exists{m : m in W | -35 eq -m and m mod 4 eq 3};   // a non-fundamental discriminant, and a point fixed by w_35
    g := 0;
    for p in [q : q in PrimesInInterval(3, 120) | D mod q ne 0] do
        s := [ p^i + 1 - ComputePointsViaTrace(Y, p, i) : i in [1..3] ];
        e := [Rationals() | 1];
        for k in [1..3] do Append(~e, (&+[ (-1)^(i-1) * e[k-i+1] * s[i] : i in [1..k] ]) / k); end for;
        // L_p(1) = 1 - e1 + e2 - e3 + p e2 - p^2 e1 + p^3
        g := GCD(g, Integers()!(1 - e[2] + e[3] - e[4] + p*e[3] - p^2*e[2] + p^3));
        if p eq 3 then assert [-e[2], e[3], -e[4]] eq [3, 9, 15]; end if;   // c1, c2, c3: T(9) has x^3 + 3x^2 + (9 - 9)x + (15 - 18)
    end for;
    assert g eq 1;
    printf " Done!\n";
end procedure;

test_hyperelliptic_4130_1();
