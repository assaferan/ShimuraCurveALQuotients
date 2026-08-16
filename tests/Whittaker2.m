// Regression test for the internal p=2 local Whittaker / representation density (Wpoly2).
//
// The value of a Borcherds form at a CM point is assembled from local Whittaker functions
// (Schofer / Kudla-Rapoport-Yang). The p=2 factor is the delicate one -- Errthum flagged an error in
// it, and it is exactly what breaks for CM points at even level (X0^D(N), 2|N), where lambda^perp is
// 2-MODULAR at 2 rather than unimodular. This test locks in that Wpoly2 agrees with Yang's explicit
// formula (Yang, "An explicit formula for local densities of quadratic forms", J. Number Theory 72
// (1998), Thm 4.1) for the off-diagonal (hyperbolic) block, in both the unimodular and 2-modular
// cases, and for integral and non-integral coset shifts.
//
// It exercises the intrinsic LocalWhittakerAtOne, which returns the (unnormalized) Whittaker
// polynomial evaluated at X=1 -- exactly the quantity kappaminus tests for vanishing.

// Independent implementation of Yang Thm 4.1 for a SINGLE off-diagonal 2-adic Jordan block
//   eps' * 2^v * [[0,1],[1,0]]   (v = 0: unimodular hyperbolic;  v = 1: 2-modular),  coset shift mu = 0.
// Yang gives  alpha(X,t,S) = 1 + R_1  with, for this block and mu = 0,
//   R_1(1) = sum_{0<k<=a+3, nu(k) in 4Z_2}  2^(min(v,k)-1) * psi(nu(k)/8),   nu(k) = m*2^(3-k), a=v_2(m).
function yang_offdiag_at_one(m, v)
    total := Rationals()!1;
    a := Valuation(m, 2);
    for k in [1..a+3] do
        nu := m * 2^(3-k);
        if Valuation(nu, 2) ge 2 then                       // char(4*Z_2)(nu)
            psi := (-1)^(Integers()!(GF(2)!(nu/4)));         // psi(nu/8) for nu in 4*Z_2
            total +:= 2^(Minimum(v, k) - 1) * psi;
        end if;
    end for;
    return total;
end function;

procedure test_Whittaker2()
    printf "Testing p=2 local Whittaker (Wpoly2) vs Yang JNT 72 (1998) Thm 4.1...";
    L := RSpaceWithBasis(IdentityMatrix(Integers(), 2));
    mu0 := Vector([Rationals() | 0, 0]);

    // (1) Off-diagonal (hyperbolic) block, mu = 0: code must equal Yang's formula over a range of m,
    //     for BOTH the unimodular (v=0) and the 2-modular (v=1, the even-level case) block.
    blocks := [ <0, Matrix(Integers(), 2, 2, [0,1,1,0])>,    // unimodular hyperbolic
                <1, Matrix(Integers(), 2, 2, [0,2,2,0])> ];   // 2-modular hyperbolic
    for bd in blocks do
        v, Qb := Explode(bd);
        for m in [1..40] do
            assert LocalWhittakerAtOne(Rationals()!m, 2, mu0, L, Qb) eq yang_offdiag_at_one(Rationals()!m, v);
        end for;
        // fractional m (2-adically a unit denominator) behaves as the integer numerator 2-adically
        for m in [1/3, 2/3, 5/7, 3/7, 10/7] do
            assert LocalWhittakerAtOne(m, 2, mu0, L, Qb) eq yang_offdiag_at_one(m, v);
        end for;
    end for;

    // (2) Hand-verified explicit values, 2-modular hyperbolic, mu = 0  (m -> W_m,2(0)):
    //     these were checked by hand against Yang Thm 4.1 and pin the exact numbers.
    Q2 := Matrix(Integers(), 2, 2, [0,2,2,0]);
    for datum in [ <1,0>, <2,1>, <3,0>, <4,2>, <5,0>, <6,1>, <7,0>, <8,3> ] do
        m, val := Explode(datum);
        assert LocalWhittakerAtOne(Rationals()!m, 2, mu0, L, Q2) eq val;
    end for;

    // (3) Non-integral (half-integral) coset shift, 2-modular hyperbolic. Here Yang's K_mu = 0, so
    //     W = char(Q(mu)+Z_2)(m): nonzero iff m is in Q(mu)+Z_2, else 0.
    //     mu=(1/2,1/2): Q(mu)=1/2, so m integer => m-Q(mu) not 2-integral => W vanishes.
    assert LocalWhittakerAtOne(Rationals()!1, 2, Vector([Rationals()|1/2,1/2]), L, Q2) eq 0;
    assert LocalWhittakerAtOne(Rationals()!3, 2, Vector([Rationals()|1/2,1/2]), L, Q2) eq 0;
    //     mu=(1/2,0): Q(mu)=0, so m integer => m in Q(mu)+Z_2 => W nonzero (= char value 1).
    assert LocalWhittakerAtOne(Rationals()!1, 2, Vector([Rationals()|1/2,0]), L, Q2) eq 1;

    printf "Done!\n";
end procedure;

test_Whittaker2();
