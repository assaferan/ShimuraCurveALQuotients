// Regression test for the local Whittaker factor at the LEVEL prime, against Kudla-Yang's closed form.
//
// Kudla-Yang, "Eisenstein series for SL(2)", Sci. China Math. 53 (2010) 2275-2316, Section 5.3 treats
// exactly the lattice this pipeline uses: eq (5.3) defines, for a quaternion algebra over Q_p with
// V the TRACE ZERO elements and Q(x) = kappa*det(x),
//
//     L_e = { (b a ; c -b) : a, b in Z_p, c in p^e Z_p },
//
// i.e. the trace-zero part of an Eichler order of level p^e.  For X0^D(N) with N squarefree and
// (N, D) = 1 this is the local lattice at p | N with e = 1.  Proposition 5.4 (B = M_2(Q_p) split,
// e > 0) then gives the local Whittaker function; it is written out in full below and compared against
// LocalWhittakerAtOne, which evaluates the Whittaker polynomial at X = 1.
//
// In the common case where the conductor exponent k_p is 0 it collapses to the memorable form
//
//     W_p(m) = 1 - 1/p          if p | d,
//     W_p(m) = 1 + chi_d(p)     if p does not divide d   (so 0 at an inert p, 2 at a split p),
//
// but that is NOT the general statement -- see the note on the conductor at ky_prop54_at_one below.
//
// WHY THIS TEST EXISTS.  Project memory long recorded that the pipeline's local factor was "a lattice
// density, not the Eisenstein coefficient", and that the ternary-with-Eichler-level case was not
// covered by the literature.  Both are wrong: Wpoly implements Kudla-Yang Theorem 4.3 (and Wpoly2
// Theorem 4.4), and Section 5.3/5.4 is exactly the ternary Eichler case.  This test pins the level
// prime against the published closed form so that identification cannot silently rot, and so that the
// m = 0 investigation cannot re-litigate whether the code is right at p | N.  (It is.)
//
// Complements tests/Whittaker2.m (p = 2 vs Yang JNT 72 (1998) Thm 4.1) and tests/M0LocalDensity.m
// (structural rep-count cross-check).  Here the comparison is against Kudla-Yang's closed form at the
// LEVEL prime specifically, on the actual ShimuraCurveLattice of real bases.

// Kudla-Yang Prop 5.4, written out independently, at X = 1 and e = 1.
//
// With c(l) = min(l/2, e - l/2) and k = k_p = ord_p(c) the CONDUCTOR exponent (4*kappa*m = d*c^2, eq
// (5.4)), Prop 5.4 reads
//     p | d      : 1 + (1-1/p) sum_{1<=l<=k} p^{c(2l)} X^{2l} - p^{c(2k+2)-1} X^{2k+2}
//     p not | d  : 1 + (1-1/p) sum_{1<=l<=k} p^{c(2l)} X^{2l} + chi_d(p) p^{c(2k+1)-1/2} X^{2k+1}
// At e = 1 and l >= 1 one has c(2l) = 1 - l, c(2k+2) = -k and c(2k+1) = 1/2 - k, so at X = 1
//     p | d      : 1 + (1-1/p) sum_{l=1..k} p^{1-l} - p^{-k-1}
//     p not | d  : 1 + (1-1/p) sum_{l=1..k} p^{1-l} + chi_d(p) p^{-k}
// (For k = 0 these collapse to 1 - 1/p and 1 + chi_d(p) respectively.)
//
// NOTE the indexing is by the CONDUCTOR, not by ord_p(m): e.g. m = 25 at p = 5 has d = -4 and c = 5,
// so k = 1 -- not the k = 0 that ord_5(25) = 2 might suggest.  Getting this wrong is an easy mistake
// and is precisely what this test guards.
// SIGN CONVENTION, pinned by this test.  The character is that of Q(sqrt(kappa m)), and for the
// Whittaker call the pipeline passes negQ = -Q, which makes the relevant field Q(sqrt(m)) -- NOT the
// Q(sqrt(-m)) used a few lines away in the same routine for the class number h/w.  The discriminating
// case is X0^6(7) at m = 1: chi_{-4}(7) = -1 would predict 0, but the code gives 2, and the trivial
// character of Q(sqrt 1) gives 2.  Verified against 26 independently dumped code values.
function ky_prop54_at_one(m, p)
    s := SquarefreeFactorization(m);
    d := (s eq 1) select 1 else Discriminant(Integers(QuadraticField(s)));
    is_sq, c := IsSquare(Rationals()!(4*m) / AbsoluteValue(d));
    error if not is_sq, Sprintf("4m/|d| is not a square for m = %o", m);
    k := Valuation(Integers()!c, p);
    chi := (d eq 1) select 1 else KroneckerSymbol(d, p);
    total := Rationals()!1;
    for l in [1..k] do
        total +:= (1 - Rationals()!1/p) * (Rationals()!p)^(1-l);
    end for;
    if d mod p eq 0 then
        total -:= (Rationals()!p)^(-k-1);
    else
        total +:= chi * (Rationals()!p)^(-k);
    end if;
    return total;
end function;

// Nonzero isotropic cosets of the FULL rank-3 lattice L^v/L that are supported only at the level
// prime N: mirrors km0places.m's coset search (Smith form of the Gram matrix) but on the full 3x3
// Q rather than the L_- restriction, and keeps only cosets whose denominators are powers of N.
function LevelIsotropicCosets(Q, N)
    Qrat := ChangeRing(Q, Rationals());
    Qz := ChangeRing(Q, Integers());
    Sm, P, R := SmithForm(Qz);
    e1 := Sm[1][1]; e2 := Sm[2][2]; e3 := Sm[3][3];
    Qi := Qrat^(-1);
    cands := [];
    for i in [0..Abs(e1)-1] do for j in [0..Abs(e2)-1] do for k in [0..Abs(e3)-1] do
        if i eq 0 and j eq 0 and k eq 0 then continue; end if;
        coef := Vector([Rationals()| i, j, k]) * ChangeRing(R^(-1), Rationals());
        cL := coef * Qi;
        cL := Vector([Rationals()| x - Floor(x) : x in Eltseq(cL)]);
        if IsZero(cL) then continue; end if;
        dens := [Denominator(x) : x in Eltseq(cL)];
        if exists{d : d in dens | d ne 1 and (d mod N ne 0 or exists{p : p in PrimeDivisors(d) | p ne N})} then
            continue;
        end if;
        q := (cL*Qrat, cL)/2;
        if not IsIntegral(q) then continue; end if;
        Append(~cands, cL);
    end for; end for; end for;
    return cands;
end function;

// Kudla-Yang Proposition 5.5 (mu in L*_e - L_e, i.e. the nonzero cosets Prop 5.4 does not cover),
// specialised to e = 1 and x1 = 0 -- the shape of the nonzero isotropic, N-only-supported cosets
// found above (one off-diagonal coordinate a unit over N, the other and the diagonal coordinate 0).
// There t_mu(m) = m - kappa*p*x2*x3 = m (since x3 = 0), so a = ord_p(t_mu(m)) is just ord_p of the
// USUAL conductor exponent k_p, and K_e = min(e + ord_p(x2), e + ord_p(x3)) = min(0, infinity) = 0.
// With K_e = 0 the sum 1 <= l <= min(k_p, K_e/2) is empty and BOTH tail conditions (2k_p+1 < 0 and
// 2k_p < 0) are vacuous for k_p >= 0, so Prop 5.5 collapses to the constant polynomial 1 -- no m
// dependence, no pole -- for every m and every conductor. This is the closed-form prediction the
// test below checks against the actual lattice.
//
// WHY THIS TEST EXISTS.  It is the mu != 0 counterpart of the mu = 0 check above (Prop 5.4).  Unlike
// that case, Prop 5.5 predicts total triviality: no m-dependence at all at the level prime, for a
// nonzero isotropic N-only-supported coset. That is a real, checkable structural claim -- and it
// falsifies the hope that inserting Prop 5.4/5.5 into a level-N analogue of Theorem 8.1 (KY section 8)
// could reproduce the b^{eta*} Eisenstein coefficients (memory: b-eisenstein-coefficients-solved),
// since those are supported exactly on N | r and vary nontrivially with N -- neither of which a
// constant local factor can produce. See PLAN.md, MAIN LINE.
procedure test_prop55_nonzero_isotropic_coset()
    printf "Testing Prop 5.5 (nonzero isotropic coset) vs LocalWhittakerAtOne...";

    bases := [ <15,2>, <6,5>, <10,3>, <21,2>, <34,11> ];
    Lfull := RSpaceWithBasis(IdentityMatrix(Integers(), 3));
    nchecked := 0;
    ncosets_total := 0;

    for b in bases do
        D, N := Explode(b);
        Ld := ShimuraCurveLattice(D, N);
        Q := ChangeRing(Ld`Q, Integers());
        negQ := -Q;

        cosets := LevelIsotropicCosets(Q, N);
        error if IsEmpty(cosets),
            Sprintf("X0^%o(%o): found no nonzero isotropic N-supported coset -- test is vacuous", D, N);
        ncosets_total +:= #cosets;

        for eta in cosets do
            error if IsZero(eta), Sprintf("X0^%o(%o): coset search returned the zero vector", D, N);
            Qeta := (eta*ChangeRing(Q, Rationals()), eta)/2;
            error if not IsIntegral(Qeta),
                Sprintf("X0^%o(%o): coset %o is not isotropic (Q = %o)", D, N, eta, Qeta);
            for m in [1..40] do
                got := LocalWhittakerAtOne(Rationals()!m, N, eta, Lfull, negQ);
                error if got ne 1,
                    Sprintf("X0^%o(%o): eta=%o m=%o gives W=%o, Kudla-Yang Prop 5.5 predicts 1",
                            D, N, eta, m, got);
                nchecked +:= 1;
            end for;
        end for;
    end for;

    printf " OK (%o checks, %o cosets over %o bases)\n", nchecked, ncosets_total, #bases;
end procedure;

procedure test_KudlaYangLocal()
    printf "Testing the level-prime local Whittaker vs Kudla-Yang Prop 5.4...";

    // (D, N) with N an odd prime not dividing D; the level prime is N.
    bases := [ <6, 5>, <6, 7>, <6, 11>, <6, 17>, <10, 3>, <10, 11>, <14, 11>, <22, 7>, <34, 5> ];
    Lfull := RSpaceWithBasis(IdentityMatrix(Integers(), 3));
    mu0 := Vector([Rationals() | 0, 0, 0]);
    nchecked := 0;

    for b in bases do
        D, N := Explode(b);
        Ld := ShimuraCurveLattice(D, N);
        negQ := -ChangeRing(Ld`Q, Integers());
        for m in [1..30] do
            got := LocalWhittakerAtOne(Rationals()!m, N, mu0, Lfull, negQ);
            expected := ky_prop54_at_one(m, N);
            error if got ne expected,
                Sprintf("X0^%o(%o): level-prime Whittaker at m=%o is %o, Kudla-Yang Prop 5.4 gives %o",
                        D, N, m, got, expected);
            nchecked +:= 1;
        end for;

        // Both branches of Prop 5.4 must actually occur, and the factor must both vanish and not
        // vanish somewhere, or the test is vacuous.
        vals := [ky_prop54_at_one(m, N) : m in [1..30]];
        error if not (exists{v : v in vals | v eq 0} and exists{v : v in vals | v ne 0}),
            Sprintf("X0^%o(%o): test is vacuous, the factor never vanishes or never does not", D, N);
    end for;

    // The distinguishing content: the factor genuinely VANISHES for some m and not others.  Pin that
    // explicitly, since the m = 0 investigation turns on which terms drop out.
    D := 6; N := 5;
    Ld := ShimuraCurveLattice(D, N);
    negQ := -ChangeRing(Ld`Q, Integers());
    vanishing := [m : m in [1..30] | LocalWhittakerAtOne(Rationals()!m, N, mu0, Lfull, negQ) eq 0];
    // With conductor exponent 0 the p-not-dividing-d branch is 1 + chi_d(p), so the factor vanishes
    // exactly at the inert primes -- for squarefree m that is the quadratic non-residues.
    expected_vanishing := [m : m in [1..30] | ky_prop54_at_one(m, N) eq 0];
    error if vanishing ne expected_vanishing,
        Sprintf("X0^6(5): level factor vanishes at %o, Kudla-Yang predicts %o",
                vanishing, expected_vanishing);
    error if IsEmpty(vanishing), "X0^6(5): vanishing check is vacuous";

    printf " OK (%o checks over %o bases)\n", nchecked, #bases;
end procedure;

test_KudlaYangLocal();
test_prop55_nonzero_isotropic_coset();
