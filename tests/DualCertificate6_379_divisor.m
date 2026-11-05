// X_0^6(379)* IS HYPERELLIPTIC, PART ONE: THE DIVISOR OF THE CERTIFICATE, THE POINTS, THE TORSION.
//
// The star C = X_0^6(379)* has genus three and carries ten rational CM points, of discriminants
// -3, -40, -43, -52, -123, -163, -168, -228, -312, -408 (each a single point: the CM points of one
// discriminant on the star form one Galois orbit, and the degree of their field of definition is 1).
// data/dual/cert_6_379.m lists the principal part c(-n), n <= 700, of a weight-1/2 form f on
// Gamma_0(4548), holomorphic away from the cusp at infinity; by Guo-Yang's Lemma 25 the Borcherds
// lift of f has divisor
//     sum_n c(-n) sum_{r : -4n/r^2 a discriminant} (1/e_d) P_d,   d = -4n/r^2,
// P_d the CM divisor of discriminant d (empty when no order of that discriminant embeds optimally
// into the Eichler order of level 379 in the quaternion algebra of discriminant 6) and e_d the
// stabiliser order of its points.  This test checks that the coefficients sum to lambda*e_d at
// d = -123, -163 (plus) and -52, -168 (minus) for one positive integer lambda, and to zero at every
// other discriminant that carries CM points, so that the lift is a function on C with divisor
// lambda*(P_-123 + P_-163 - P_-52 - P_-168).  Its class is therefore torsion; Jac(C)(Q) is
// torsion-free (the gcd of #Jac(C)(F_p) over good primes is 1, from the trace formula and no model),
// so P_-123 + P_-163 ~ P_-52 + P_-168: two distinct effective divisors of degree two in one class,
// i.e. a g^1_2, and C is hyperelliptic.  The existence of f (the pairing with the weight-3/2 cusp
// forms of level 4548) is tests/_offline/DualCertificate6_379.m.

// the certificate is read at top level: evaluating a data file inside a procedure crashes Magma here
cert := eval (Read("data/dual/cert_6_379.m") cat " return cert;");

procedure test_dual_certificate_6_379_divisor(cert)
    printf "X_0^6(379)* is hyperelliptic: the certificate's divisor is lambda(P_-123 + P_-163 - P_-52 - P_-168) and Jac(Q) is torsion-free...";
    D := 6; N := 379;
    W := {m : m in Divisors(D*N) | GCD(m, (D*N) div m) eq 1};
    Xstar := CreateShimuraQuot(D, N, W);
    Xstar`g := GenusShimuraCurveQuotient(D, N, W); Xstar`CurveID := 0;
    assert Xstar`g eq 3;
    curves := GetQuotientsAndGenera([Xstar]);
    _ := exists(star){c : c in curves | IsStarCurve(c)};

    // (1) the four points are rational points of the star, one per discriminant
    rat := RationalandQuadraticCMPoints(star : bd := 4, coprime_to_level := false);
    assert {p[1] : p in rat} eq {-3, -40, -43, -52, -123, -163, -168, -228, -312, -408};
    assert forall{d : d in [-123, -163, -52, -168] | DegreeOfFieldOfDefinitionOfCMPoint(star, d) eq 1};

    // (2) the divisor of the lift, from the principal part
    B := #cert;
    function HasCM(d)
        OK := MaximalOrder(QuadraticField(d)); _, f := IsSquare(d div Discriminant(OK));
        return NumberOfOptimalEmbeddings(sub<OK | f>, D, N) gt 0;
    end function;
    function Stab(d)
        e := d eq -3 select 3 else (d eq -4 select 2 else 1);
        if exists{m : m in Divisors(D*N) | m gt 1 and GCD(m, (D*N) div m) eq 1 and (d eq -4*m or (m mod 4 eq 3 and d eq -m))} then e *:= 2; end if;
        return e;
    end function;
    target := AssociativeArray(); target[-123] := 1; target[-163] := 1; target[-52] := -1; target[-168] := -1;
    sums := AssociativeArray();
    for n in [1..B] do
        if cert[n] eq 0 then continue; end if;
        for r in [1..Isqrt(4*n)] do
            if (4*n) mod (r^2) ne 0 then continue; end if;
            d := -4*n div (r^2);
            if (d mod 4) notin {0, 1} then continue; end if;
            if not HasCM(d) then continue; end if;
            sums[d] := (IsDefined(sums, d) select sums[d] else 0) + cert[n];
        end for;
    end for;
    // the divisor is lambda times the target for one positive integer lambda (the vector is integral
    // and primitive; a multiple of the divisor is all the argument needs)
    assert forall{ d : d in Keys(target) | IsDefined(sums, d) };
    lambda := -sums[-52];
    assert lambda gt 0 and IsIntegral(lambda);
    bad := [ d : d in Keys(sums) | sums[d] ne (IsDefined(target, d) select lambda*Stab(d)*target[d] else 0) ];
    assert IsEmpty(bad);
    assert forall{d : d in Keys(target) | Stab(d) eq 1};
    assert GCD([Integers() | c : c in cert]) eq 1;

    // (3) Jac(C)(Q) is torsion-free: the torsion injects into Jac(C)(F_p) for odd good p, and
    // #Jac(C)(F_p) = P(1) for the characteristic polynomial P of Frobenius, read off the point
    // counts over F_p, F_{p^2}, F_{p^3} from the trace formula (genus three: Newton's identities
    // and the functional equation)
    t := 0;
    for p in [q : q in PrimesInInterval(5, 60) | (D*N) mod q ne 0] do
        s := [ p^i + 1 - ComputePointsViaTrace(Xstar, p, i) : i in [1..3] ];
        e1 := s[1]; e2 := (e1*s[1] - s[2]) div 2; e3 := (e2*s[1] - e1*s[2] + s[3]) div 3;
        assert 2*e2 eq e1*s[1] - s[2] and 3*e3 eq e2*s[1] - e1*s[2] + s[3];
        nJ := 1 - e1 + e2 - e3 + p*e2 - p^2*e1 + p^3;
        t := GCD(t, nJ);
    end for;
    assert t eq 1;
    printf " [%o CM discriminants met, %o nonzero coefficients, pole order %o, lambda = %o, torsion bound %o] Done!\n",
        #Keys(sums), #[n : n in [1..B] | cert[n] ne 0], Maximum([n : n in [1..B] | cert[n] ne 0]), lambda, t;
end procedure;

test_dual_certificate_6_379_divisor(cert);
