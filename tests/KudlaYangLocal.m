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
