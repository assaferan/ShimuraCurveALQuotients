// THE DUAL CERTIFICATE OF X_0^2002(1)/<w_2,w_13,w_77>, PART ONE: ITS DIVISOR.
//
// data/dual/cert_2002_1.m lists the principal part c(-n), n <= 2200, of a weight-1/2 form f on
// Gamma_0(4004).  By Guo-Yang's Lemma 25 the Borcherds lift of f has divisor
//     sum_n c(-n) sum_{r : -4n/r^2 a discriminant} (1/e_d) P_d,   d = -4n/r^2,
// on the star, P_d the CM divisor of discriminant d (empty when no order of that discriminant embeds
// optimally into the maximal order of discriminant 2002) and e_d the stabiliser order of its points
// (2 at an Atkin-Lehner fixed point).  This test checks that the coefficients sum to e_d at d = -11,
// -91 (plus one each) and -148, -232 (minus one each) and to zero at every other discriminant that
// carries CM points, so that the divisor is P_-11 + P_-91 - P_-148 - P_-232.  The existence of f
// (the pairing with the weight-3/2 cusp forms) is tests/_offline/DualCertificate2002_1.m.

// the certificate is read at top level: evaluating a data file inside a procedure crashes Magma here
cert := eval (Read("data/dual/cert_2002_1.m") cat " return cert;");

procedure test_dual_certificate_2002_1_divisor(cert)
    printf "the dual certificate of X_0^2002(1)/<w_2,w_13,w_77> has divisor P_-11 + P_-91 - P_-148 - P_-232...";
    D := 2002; N := 1;
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
    target := AssociativeArray(); target[-11] := 1; target[-91] := 1; target[-148] := -1; target[-232] := -1;
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
    // and primitive; a multiple of the divisor is all the counting argument needs)
    assert forall{ d : d in Keys(target) | IsDefined(sums, d) };
    lambda := -sums[-148];
    assert lambda gt 0 and IsIntegral(lambda);
    bad := [ d : d in Keys(sums) | sums[d] ne (IsDefined(target, d) select lambda*Stab(d)*target[d] else 0) ];
    assert IsEmpty(bad);
    assert Stab(-11) eq 2 and Stab(-91) eq 2 and Stab(-148) eq 1 and Stab(-232) eq 1;
    assert GCD([Integers() | c : c in cert]) eq 1;
    printf " [%o CM discriminants met, %o nonzero coefficients, pole order %o, lambda = %o] Done!\n", #Keys(sums), #[n : n in [1..B] | cert[n] ne 0], Maximum([n : n in [1..B] | cert[n] ne 0]), lambda;
end procedure;

test_dual_certificate_2002_1_divisor(cert);
