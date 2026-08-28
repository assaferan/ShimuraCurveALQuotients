// kappa-minus, step 1: an INDEPENDENT local density, calibrated where the code
// is already trusted.
//
// tests/Whittaker2.m cross-checks LocalWhittakerAtOne against a hand-written
// Yang Thm 4.1 -- but only at p = 2, only hyperbolic blocks, only v = 0,1.  At a
// firing discriminant with N > 1 the binary L_- is N-SCALED at the level prime,
// and for ODD N (3,5,7,11,... -- most N>1 bases) nothing independent has ever
// been compared against.  That is the gap the preprint flags.
//
// The independent object is the naive local representation density
//     alpha_p(m) = lim_k p^{-k} #{ x in (Z/p^k)^2 : Q(x) = m mod p^k }
// (rank 2, so the count grows like p^k).  Convergence in k is checked, not
// assumed.  Normalisation between alpha and W(1) is CALIBRATED at p = 2 on the
// blocks Whittaker2.m already validates, then carried to the untested cases; and
// since any leftover constant is m-independent, the comparison is done on RATIOS
// W(m)/W(m_ref) against alpha(m)/alpha(m_ref), which no normalisation can fake.
AttachSpec("ShimuraQuotients.spec");

// brute-force count, with an explicit convergence check in k
alpha := function(p, Qm, m, kmax)
    prev := Rationals()!0; got := false; val := Rationals()!0;
    for k in [2..kmax] do
        pk := p^k; cnt := 0;
        mm := Integers()!(m mod pk);
        for x1 in [0..pk-1] do
            for x2 in [0..pk-1] do
                q := (Qm[1][1]*x1^2 + 2*Qm[1][2]*x1*x2 + Qm[2][2]*x2^2);
                if (q - mm) mod pk eq 0 then cnt +:= 1; end if;
            end for;
        end for;
        cur := Rationals()!cnt / pk;
        if k ge 3 and cur eq prev then got := true; val := cur; break; end if;
        prev := cur;
    end for;
    return val, got;
end function;

L := RSpaceWithBasis(IdentityMatrix(Integers(), 2));
mu0 := Vector([Rationals() | 0, 0]);

// Q(x) convention probe: the code takes a Gram matrix; try x^T Q x and half of it
printf "=== calibration at p = 2, hyperbolic (the blocks Whittaker2.m validates) ===\n";
for vv in [0, 1] do
    Qb := Matrix(Integers(), 2, 2, [0, 2^vv, 2^vv, 0]);
    printf "  v = %o,  Q = %o\n", vv, Eltseq(Qb);
    for m in [1, 2, 3, 4, 5, 6, 8, 12] do
        w := LocalWhittakerAtOne(Rationals()!m, 2, mu0, L, Qb);
        a1, ok1 := alpha(2, Qb, m, 9);            // Q(x) = x^T Q x
        printf "    m=%-3o  W(1) = %-8o  alpha[xQx] = %-10o %o  ratio %o\n",
            m, w, a1, ok1 select "" else "(NOT CONVERGED)",
            a1 ne 0 select Sprintf("%o", w/a1) else "-";
    end for;
end for;
printf "DONE\n";
quit;
