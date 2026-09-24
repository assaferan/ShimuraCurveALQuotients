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
//
// ---------------------------------------------------------------------------------------------
// WHICH OBJECTS THIS COVERS, AND WHY THESE ONES  (parts (4) and (5), added 2026-09-15)
//
// Wpoly2 is on the production path: Wpoly_scaled dispatches to it at p = 2, feeding get_Wpolys ->
// get_kappa_minus_squared -> the Schofer CM values. Sweeping 9 bases x 10 discriminants and
// 2-adically diagonalising the Gram of lambda^perp, exactly FIVE 2-adic Jordan shapes occur as
// production input -- and parts (1)-(3) covered only two of them:
//
//     even0H           unimodular hyperbolic  [[0,1],[1,0]]            <- parts (1)-(3)
//     even1H           2-modular hyperbolic   [[0,2],[2,0]]            <- parts (1)-(3)
//     d1d1             diagonal, 2-valuations 1,1  (v_2(det) = 2)      <- parts (4),(5)
//     d1d2             diagonal, 2-valuations 1,2  (v_2(det) = 3)      <- parts (4),(5)
//     even1A           2-modular NON-hyperbolic  2*[[2,1],[1,2]]       <- parts (4),(5)
//
// None of the three new ones is exotic. d1d1/d1d2: 6_1 at d = -4 and -24, 14_1 at -4 and -8, 10_1 at
// -8 and -20, 34_1 at -20 and -24, 38_1 at -4/-20/-24, 26_1 at -8/-20/-24, 10_3 at -8/-20, 6_5 at
// -4/-24, 21_1 at -4. even1A: it is what the ODD discriminants give -- 6_1 at -3 and -19, 14_1 at
// -11, 10_1 at -3, 34_1 at -3/-11, 38_1 at -11/-19, 26_1 at -11/-19, 6_5 at -19, 10_3 at -3. All
// twelve of those, over eight different bases, land in the SAME 2-adic class as 2*[[2,1],[1,2]] and
// return the same values.
//
// Until 2026-09-15 nothing tested any of the three: parts (1)-(3) only ever hand Wpoly2 an even
// HYPERBOLIC block. even1A is even too but not hyperbolic, and Wpoly2 routes it through the
// x^2 + xy + y^2 blocks instead -- it returns values (3, and 0 at m = 4) that neither hyperbolic
// block ever produces, so it really is a separate branch and not a relabelling.
//
// HOW THE EXPECTED VALUES IN (4) AND (5) WERE OBTAINED. Not from Yang's formula and not from this
// library -- from a brute-force representation-density count, which is independent of both:
//
//     alpha_2(m, L) = lim_{k->oo} 2^(-k(n-1)) * #{ x in (mu + Z/2^k)^n : Q(x) = m mod 2^k },
//
// with this repo's convention Q(x) = 1/2 * x G x^T, i.e. counting x G x^T = 2m mod 2^(k+1). Run
// offline (tests/_offline/Whittaker2Oracle.m regenerates the whole table); the committed numbers are
// the CONVERGED limits. Calibration first, per this repo's standing rule: the same counter
// reproduces LocalWhittakerAtOne EXACTLY, ratio 1, on 32/32 of parts (1)-(3)'s own configurations:
// m in [1,2,3,4,5,6,8,12] on each hyperbolic block at mu = 0, and the same eight m on each of part
// (3)'s two non-integral cosets. If it does not, the counter is wrong -- fix it before going on.
//
// *** THE CONVERGENCE TRAP -- DO NOT "OPTIMISE" THE ORACLE WITH A REPEATS-BASED STOPPING RULE. ***
// The approximants PLATEAU AT A WRONG VALUE before converging, and the plateau is long. Measured at
// 6_1, d = -24 (shape d1d2), k = 1..14:
//
//     m = 8    1 1 1 1 1 2 2 2 2 2 2 2 2 2      (true value 2, reached only at k = 6)
//     m = 32   1 1 1 1 1 1 1 2 2 2 2 2 2 2      (true value 2, reached only at k = 8)
//
// Accepting "two successive repeats" as convergence returns 1 for both and produces false
// mismatches, concentrated at the most 2-divisible m. There is no repeats rule here: the table was
// taken at a FIXED k = 14 and re-confirmed at k = 16. Across all 300 mu = 0 comparisons the LAST k
// at which any approximant still moved was k = 8; across the 666 nonzero-coset comparisons it was
// k = 2. So k = 14 carries six spare rungs at worst, and k = 16 carries eight.
// ---------------------------------------------------------------------------------------------

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

// The 2-adic Jordan shape of an integral binary Gram G, as a label.
// Let s = min over ALL entries of v_2, the 2-adic SCALE. The lattice is of EVEN type at 2 exactly
// when no DIAGONAL entry attains s; otherwise it is of ODD type, i.e. 2-adically diagonalisable as
// <2^i u> + <2^j u'> with i = s and i + j = v_2(det G) -- so s = 1 with v_2(det) = 2 is "d1d1" and
// s = 1 with v_2(det) = 3 is "d1d2".
// There are exactly TWO even binary lattices over Z_2 up to scaling: H = [[0,1],[1,0]] (det -1) and
// A = [[2,1],[1,2]] (det 3). Scaling by 2^s multiplies the determinant by 4^s, so det/4^s is a unit
// and its class mod 8 separates them: 7 is H, 3 is A. They are NOT interchangeable -- Wpoly2 sends
// them down different branches of Yang's formula (the mu' product blocks vs the mu'' x^2+xy+y^2
// blocks) and they give different values, so "even1" alone would be a WRONG-OBJECT label.
// This is the which-object guard on the committed Gram matrices below: a typo in one of them almost
// certainly changes its shape, and then this fails before any value is compared.
function two_adic_shape(G)
    v := func< x | x eq 0 select 10^6 else Valuation(x, 2) >;
    s  := Minimum([v(G[1,1]), v(G[1,2]), v(G[2,2])]);
    dg := Minimum([v(G[1,1]), v(G[2,2])]);
    if dg gt s then
        u := Determinant(G) div 4^s;
        assert u mod 8 in {3, 7};
        return Sprintf("even%o%o", s, (u mod 8 eq 7) select "H" else "A");
    end if;
    return Sprintf("d%od%o", s, Valuation(Determinant(G), 2) - s);
end function;

// The 2-primary part of L^dual/L for an integral binary Gram G, as coset representatives in the
// standard basis. An element of the 2-primary part has order dividing 2^t with t = v_2(det G), so
// its coordinates have denominators dividing 2^t; conversely such a vector lies in the 2-primary
// part as soon as it lies in L^dual. #(L^dual/L)_2 = 2^t, which the caller asserts.
function two_primary_cosets(G)
    t  := Valuation(Determinant(G), 2);
    pt := 2^t;
    GQ := ChangeRing(G, Rationals());
    reps := [];
    for a in [0..pt-1], b in [0..pt-1] do
        w := Vector([Rationals() | a/pt, b/pt]);
        if &and[IsIntegral(x) : x in Eltseq(w*GQ)] then
            Append(~reps, w);
        end if;
    end for;
    return reps;
end function;

// Real lambda^perp Gram matrices from production, one row per (base, discriminant), harvested with
//     Ld := ShimuraCurveLattice(D,N);  Q := ChangeRing(Ld`Qinv^(-1), Integers());
//     lam := ElementOfNorm(Q, -d, Ld`O, Ld`basis_L);  Lm := Kernel(Transpose(Matrix(lam*Q)));
//     G := Matrix(Basis(Lm))*Q*Transpose(Matrix(Basis(Lm)));
// They are LITERALS, not re-derived, on purpose: ElementOfNorm consumes randomness, so which lambda
// of a given norm it returns shifts when anything upstream of it changes, and two runs of the
// harvesting script above already produced two different (2-adically equivalent) Gram matrices at
// 10_3, d = -3. What is pinned here is the object -- a real 2-adic shape that production feeds
// Wpoly2 -- and two_adic_shape below re-checks that each row still IS that shape.
//
// Row: < base, d, shape, G as [a,b,b,c],
//        [ W_{m,2}(0) for m in ms0 ],
//        [ <mu, [ W_{m,2}(mu) for m in Q(mu) + {0,1,2,3,4,8} ]> for each NONZERO 2-primary coset ] >
production_shapes := [*
  < "6_1", -4, "d1d1", [-6,-114,-114,-2172],
    [2,2,0,2,2,0,0,2,0,2,0,2],
    [ <[Rationals()|1/2,1/2], [1,1,1,1,1,1]>, <[Rationals()|0,1/2], [2,0,2,0,2,2]>, <[Rationals()|1/2,0], [1,1,1,1,1,1]> ] >,
  < "14_1", -4, "d1d1", [-1414,-2968,-2968,-6230],
    [2,2,0,2,2,0,0,2,0,2,0,2],
    [ <[Rationals()|0,1/2], [1,1,1,1,1,1]>, <[Rationals()|1/2,0], [1,1,1,1,1,1]>, <[Rationals()|1/2,1/2], [2,0,2,0,2,2]> ] >,
  < "10_1", -20, "d1d1", [-326,370,370,-420],
    [2,0,0,2,2,2,0,0,0,2,2,0],
    [ <[Rationals()|1/2,1/2], [1,1,1,1,1,1]>, <[Rationals()|0,1/2], [0,2,0,2,0,0]>, <[Rationals()|1/2,0], [1,1,1,1,1,1]> ] >,
  < "21_1", -4, "d1d1", [-42,0,0,-42],
    [0,0,2,0,0,2,2,0,2,0,2,0],
    [ <[Rationals()|0,1/2], [1,1,1,1,1,1]>, <[Rationals()|1/2,0], [1,1,1,1,1,1]>, <[Rationals()|1/2,1/2], [0,2,0,2,0,0]> ] >,
  < "26_1", -20, "d1d1", [-78,416,416,-2262],
    [2,0,0,2,2,2,0,0,0,2,2,0],
    [ <[Rationals()|0,1/2], [1,1,1,1,1,1]>, <[Rationals()|1/2,1/2], [0,2,0,2,0,0]>, <[Rationals()|1/2,0], [1,1,1,1,1,1]> ] >,
  < "34_1", -20, "d1d1", [-33558,10574,10574,-3332],
    [2,0,0,2,2,2,0,0,0,2,2,0],
    [ <[Rationals()|1/2,1/2], [1,1,1,1,1,1]>, <[Rationals()|0,1/2], [0,2,0,2,0,0]>, <[Rationals()|1/2,0], [1,1,1,1,1,1]> ] >,
  < "38_1", -4, "d1d1", [-52630,-92454,-92454,-162412],
    [2,2,0,2,2,0,0,2,0,2,0,2],
    [ <[Rationals()|1/2,1/2], [1,1,1,1,1,1]>, <[Rationals()|0,1/2], [2,0,2,0,2,2]>, <[Rationals()|1/2,0], [1,1,1,1,1,1]> ] >,
  < "38_1", -20, "d1d1", [-342,-2470,-2470,-17860],
    [2,0,0,2,2,2,0,0,0,2,2,0],
    [ <[Rationals()|1/2,1/2], [1,1,1,1,1,1]>, <[Rationals()|0,1/2], [0,2,0,2,0,0]>, <[Rationals()|1/2,0], [1,1,1,1,1,1]> ] >,
  < "10_3", -20, "d1d1", [-486,-2310,-2310,-10980],
    [2,0,0,2,2,2,0,0,0,2,2,0],
    [ <[Rationals()|1/2,1/2], [1,1,1,1,1,1]>, <[Rationals()|0,1/2], [0,2,0,2,0,0]>, <[Rationals()|1/2,0], [1,1,1,1,1,1]> ] >,
  < "6_5", -4, "d1d1", [-870,-510,-510,-300],
    [2,2,0,2,2,0,0,2,0,2,0,2],
    [ <[Rationals()|1/2,1/2], [1,1,1,1,1,1]>, <[Rationals()|0,1/2], [2,0,2,0,2,2]>, <[Rationals()|1/2,0], [1,1,1,1,1,1]> ] >,
  < "6_1", -24, "d1d2", [-22,-52,-52,-124],
    [0,2,2,0,2,0,0,2,2,0,0,2],
    [ <[Rationals()|1/2,1/4], [1,1,1,1,1,1]>, <[Rationals()|0,1/2], [2,0,0,2,2,2]>, <[Rationals()|1/2,3/4], [1,1,1,1,1,1]>, <[Rationals()|0,1/4], [1,1,1,1,1,1]>, <[Rationals()|1/2,1/2], [2,0,2,0,2,2]>, <[Rationals()|0,3/4], [1,1,1,1,1,1]>, <[Rationals()|1/2,0], [0,2,0,2,0,0]> ] >,
  < "14_1", -8, "d1d2", [-714,-1568,-1568,-3444],
    [2,2,2,2,0,2,0,2,2,2,2,2],
    [ <[Rationals()|0,1/4], [1,1,1,1,1,1]>, <[Rationals()|0,1/2], [2,2,0,0,2,2]>, <[Rationals()|0,3/4], [1,1,1,1,1,1]>, <[Rationals()|1/2,0], [2,0,2,0,2,2]>, <[Rationals()|1/2,1/4], [1,1,1,1,1,1]>, <[Rationals()|1/2,1/2], [2,0,2,0,2,2]>, <[Rationals()|1/2,3/4], [1,1,1,1,1,1]> ] >,
  < "10_1", -8, "d1d2", [-5930,2600,2600,-1140],
    [2,2,2,2,0,2,0,2,2,2,2,2],
    [ <[Rationals()|0,3/4], [1,1,1,1,1,1]>, <[Rationals()|0,1/2], [2,2,0,0,2,2]>, <[Rationals()|0,1/4], [1,1,1,1,1,1]>, <[Rationals()|1/2,0], [2,0,2,0,2,2]>, <[Rationals()|1/2,3/4], [1,1,1,1,1,1]>, <[Rationals()|1/2,1/2], [2,0,2,0,2,2]>, <[Rationals()|1/2,1/4], [1,1,1,1,1,1]> ] >,
  < "26_1", -8, "d1d2", [-26,234,234,-2158],
    [2,2,2,2,0,2,0,2,2,2,2,2],
    [ <[Rationals()|3/4,3/4], [1,1,1,1,1,1]>, <[Rationals()|1/2,1/2], [2,2,0,0,2,2]>, <[Rationals()|1/4,1/4], [1,1,1,1,1,1]>, <[Rationals()|1/4,3/4], [1,1,1,1,1,1]>, <[Rationals()|0,1/2], [2,0,2,0,2,2]>, <[Rationals()|3/4,1/4], [1,1,1,1,1,1]>, <[Rationals()|1/2,0], [2,0,2,0,2,2]> ] >,
  < "26_1", -24, "d1d2", [-650,988,988,-1508],
    [0,2,2,0,2,0,0,2,2,0,0,2],
    [ <[Rationals()|1/2,3/4], [1,1,1,1,1,1]>, <[Rationals()|0,1/2], [2,0,0,2,2,2]>, <[Rationals()|1/2,1/4], [1,1,1,1,1,1]>, <[Rationals()|0,3/4], [1,1,1,1,1,1]>, <[Rationals()|1/2,1/2], [0,2,0,2,0,0]>, <[Rationals()|0,1/4], [1,1,1,1,1,1]>, <[Rationals()|1/2,0], [2,0,2,0,2,2]> ] >,
  < "34_1", -24, "d1d2", [-65450,2108,2108,-68],
    [0,2,2,0,2,0,0,2,2,0,0,2],
    [ <[Rationals()|1/2,3/4], [1,1,1,1,1,1]>, <[Rationals()|0,1/2], [2,0,0,2,2,2]>, <[Rationals()|1/2,1/4], [1,1,1,1,1,1]>, <[Rationals()|1/2,0], [2,0,2,0,2,2]>, <[Rationals()|0,3/4], [1,1,1,1,1,1]>, <[Rationals()|1/2,1/2], [0,2,0,2,0,0]>, <[Rationals()|0,1/4], [1,1,1,1,1,1]> ] >,
  < "38_1", -24, "d1d2", [-266,-836,-836,-2660],
    [0,2,2,0,2,0,0,2,2,0,0,2],
    [ <[Rationals()|1/2,1/4], [1,1,1,1,1,1]>, <[Rationals()|0,1/2], [2,0,0,2,2,2]>, <[Rationals()|1/2,3/4], [1,1,1,1,1,1]>, <[Rationals()|0,1/4], [1,1,1,1,1,1]>, <[Rationals()|1/2,1/2], [0,2,0,2,0,0]>, <[Rationals()|0,3/4], [1,1,1,1,1,1]>, <[Rationals()|1/2,0], [2,0,2,0,2,2]> ] >,
  < "10_3", -8, "d1d2", [-9690,-15000,-15000,-23220],
    [2,2,2,2,0,2,0,2,2,2,2,2],
    [ <[Rationals()|0,1/4], [1,1,1,1,1,1]>, <[Rationals()|0,1/2], [2,2,0,0,2,2]>, <[Rationals()|0,3/4], [1,1,1,1,1,1]>, <[Rationals()|1/2,0], [2,0,2,0,2,2]>, <[Rationals()|1/2,1/4], [1,1,1,1,1,1]>, <[Rationals()|1/2,1/2], [2,0,2,0,2,2]>, <[Rationals()|1/2,3/4], [1,1,1,1,1,1]> ] >,
  < "6_5", -24, "d1d2", [-10,-220,-220,-4900],
    [0,2,2,0,2,0,0,2,2,0,0,2],
    [ <[Rationals()|1/2,3/4], [1,1,1,1,1,1]>, <[Rationals()|0,1/2], [2,0,0,2,2,2]>, <[Rationals()|1/2,1/4], [1,1,1,1,1,1]>, <[Rationals()|0,3/4], [1,1,1,1,1,1]>, <[Rationals()|1/2,1/2], [0,2,0,2,0,0]>, <[Rationals()|0,1/4], [1,1,1,1,1,1]>, <[Rationals()|1/2,0], [2,0,2,0,2,2]> ] >,
  // even1A: the 2-modular NON-hyperbolic block, 2*[[2,1],[1,2]] up to Z_2-equivalence. Six real
  // Gram matrices from six different bases, all in this one class. The mu = 0 row is the
  // discriminating one -- the value 3 appears nowhere else in this file. The nonzero cosets are a
  // WEAK check here: Yang's K_mu vanishes on all three, so they only test the support rule.
  < "6_1", -3, "even1A", [-964,-62,-62,-4],
    [0,3,0,0,0,3,0,3,0,0,3,3],
    [ <[Rationals()|1/2,0], [1,1,1,1,1,1]>, <[Rationals()|0,1/2], [1,1,1,1,1,1]>, <[Rationals()|1/2,1/2], [1,1,1,1,1,1]> ] >,
  < "14_1", -11, "even1A", [-252,-1694,-1694,-11396],
    [0,3,0,0,0,3,0,3,0,0,3,3],
    [ <[Rationals()|1/2,0], [1,1,1,1,1,1]>, <[Rationals()|0,1/2], [1,1,1,1,1,1]>, <[Rationals()|1/2,1/2], [1,1,1,1,1,1]> ] >,
  < "10_1", -3, "even1A", [-9780,1170,1170,-140],
    [0,3,0,0,0,3,0,3,0,0,3,3],
    [ <[Rationals()|1/2,0], [1,1,1,1,1,1]>, <[Rationals()|0,1/2], [1,1,1,1,1,1]>, <[Rationals()|1/2,1/2], [1,1,1,1,1,1]> ] >,
  < "26_1", -19, "even1A", [-52,1430,1430,-39572],
    [0,3,0,0,0,3,0,3,0,0,3,3],
    [ <[Rationals()|1/2,0], [1,1,1,1,1,1]>, <[Rationals()|0,1/2], [1,1,1,1,1,1]>, <[Rationals()|1/2,1/2], [1,1,1,1,1,1]> ] >,
  < "34_1", -11, "even1A", [-139740,5338,5338,-204],
    [0,3,0,0,0,3,0,3,0,0,3,3],
    [ <[Rationals()|1/2,0], [1,1,1,1,1,1]>, <[Rationals()|0,1/2], [1,1,1,1,1,1]>, <[Rationals()|1/2,1/2], [1,1,1,1,1,1]> ] >,
  < "38_1", -19, "even1A", [-20,-2,-2,-4],
    [0,3,0,0,0,3,0,3,0,0,3,3],
    [ <[Rationals()|1/2,0], [1,1,1,1,1,1]>, <[Rationals()|0,1/2], [1,1,1,1,1,1]>, <[Rationals()|1/2,1/2], [1,1,1,1,1,1]> ] >

*];

procedure test_Whittaker2()
    printf "Testing p=2 local Whittaker (Wpoly2) vs Yang JNT 72 Thm 4.1 and vs brute-force densities...";
    L := RSpaceWithBasis(IdentityMatrix(Integers(), 2));
    mu0 := Vector([Rationals() | 0, 0]);
    ncmp := 0;          // every value comparison made, counted; asserted at the end

    // (1) Off-diagonal (hyperbolic) block, mu = 0: code must equal Yang's formula over a range of m,
    //     for BOTH the unimodular (v=0) and the 2-modular (v=1, the even-level case) block.
    blocks := [ <0, Matrix(Integers(), 2, 2, [0,1,1,0])>,    // unimodular hyperbolic
                <1, Matrix(Integers(), 2, 2, [0,2,2,0])> ];   // 2-modular hyperbolic
    for bd in blocks do
        v, Qb := Explode(bd);
        assert two_adic_shape(Qb) eq Sprintf("even%oH", v);
        for m in [1..40] do
            assert LocalWhittakerAtOne(Rationals()!m, 2, mu0, L, Qb) eq yang_offdiag_at_one(Rationals()!m, v);
            ncmp +:= 1;
        end for;
        // fractional m (2-adically a unit denominator) behaves as the integer numerator 2-adically
        for m in [1/3, 2/3, 5/7, 3/7, 10/7] do
            assert LocalWhittakerAtOne(m, 2, mu0, L, Qb) eq yang_offdiag_at_one(m, v);
            ncmp +:= 1;
        end for;
    end for;

    // (2) Hand-verified explicit values, 2-modular hyperbolic, mu = 0  (m -> W_m,2(0)):
    //     these were checked by hand against Yang Thm 4.1 and pin the exact numbers.
    Q2 := Matrix(Integers(), 2, 2, [0,2,2,0]);
    for datum in [ <1,0>, <2,1>, <3,0>, <4,2>, <5,0>, <6,1>, <7,0>, <8,3> ] do
        m, val := Explode(datum);
        assert LocalWhittakerAtOne(Rationals()!m, 2, mu0, L, Q2) eq val;
        ncmp +:= 1;
    end for;

    // (3) Non-integral (half-integral) coset shift, 2-modular hyperbolic. Here Yang's K_mu = 0, so
    //     W = char(Q(mu)+Z_2)(m): nonzero iff m is in Q(mu)+Z_2, else 0.
    //     mu=(1/2,1/2): Q(mu)=1/2, so m integer => m-Q(mu) not 2-integral => W vanishes.
    assert LocalWhittakerAtOne(Rationals()!1, 2, Vector([Rationals()|1/2,1/2]), L, Q2) eq 0;
    assert LocalWhittakerAtOne(Rationals()!3, 2, Vector([Rationals()|1/2,1/2]), L, Q2) eq 0;
    //     mu=(1/2,0): Q(mu)=0, so m integer => m in Q(mu)+Z_2 => W nonzero (= char value 1).
    assert LocalWhittakerAtOne(Rationals()!1, 2, Vector([Rationals()|1/2,0]), L, Q2) eq 1;
    ncmp +:= 3;

    // ------------------------------------------------------------------------------------------
    // (4) and (5): the three production shapes parts (1)-(3) miss -- d1d1, d1d2 and even1A -- on the
    //     real lambda^perp Gram matrices listed above. Expected values are the converged brute-force
    //     densities; see the file header for how they were obtained and for the convergence trap.
    ms0 := [1,2,3,4,5,6,7,8,12,16,24,32];
    nshape := AssociativeArray();
    nshape["d1d1"] := 0; nshape["d1d2"] := 0; nshape["even1A"] := 0;
    for row in production_shapes do
        base, d, shape, gent, w0, wmus := Explode(row);
        G := MatrixAlgebra(Integers(), 2) ! gent;       // AlgMatElt: LocalWhittakerAtOne needs this

        // WHICH OBJECT: this Gram really is the claimed 2-adic shape, and nothing else.
        assert two_adic_shape(G) eq shape;
        assert shape in {"d1d1", "d1d2", "even1A"};
        nshape[shape] +:= 1;

        // (4) mu = 0.
        assert #w0 eq #ms0;
        for i in [1..#ms0] do
            assert LocalWhittakerAtOne(Rationals()!ms0[i], 2, mu0, L, G) eq w0[i];
            ncmp +:= 1;
        end for;

        // (5) every NONZERO coset of the 2-primary part of L^dual/L. For these shapes Q(mu) is never
        //     2-integral at such a mu, so the interesting m are the ones in Q(mu) + Z -- at integral
        //     m the answer is 0 by support alone and tests nothing about Yang's formula.
        cosets := two_primary_cosets(G);
        assert #cosets eq 2^Valuation(Determinant(G), 2);          // = #(L^dual/L)_2
        assert {c : c in cosets} eq {mu0} join {Vector(t[1]) : t in wmus};   // all of them, none twice
        assert #wmus eq #cosets - 1;
        for t in wmus do
            mu, expected := Explode(t);
            muv := Vector(mu);
            Qmu := (Matrix(muv)*ChangeRing(G, Rationals())*Transpose(Matrix(muv)))[1,1] / 2;
            frac := Qmu - Floor(Qmu);
            assert frac ne 0;                                       // else the m list below is wrong
            ms := [frac + j : j in [0,1,2,3,4,8]];
            assert #expected eq #ms;
            for i in [1..#ms] do
                assert LocalWhittakerAtOne(ms[i], 2, muv, L, G) eq expected[i];
                ncmp +:= 1;
            end for;
        end for;
    end for;

    // COUNT THE CHECKS. A silently emptied loop -- an unpopulated table, a shape filter that matches
    // nothing, a mu enumeration that returns only the zero coset -- passes every assertion above and
    // is caught only here.
    assert nshape["d1d1"] eq 10 and nshape["d1d2"] eq 9 and nshape["even1A"] eq 6;
    assert ncmp eq 1067;    // 90 + 8 + 3 (parts 1-3) + 300 (part 4) + 666 (part 5)
    printf " Done!  %o comparisons.\n", ncmp;
end procedure;

test_Whittaker2();
