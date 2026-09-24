import "tests/BorcherdsProducts.m" : test_AllEquationsAboveCoversSingleCurve;

// tests/X0_21_2.m -- RE-DERIVATION test for X_0^21(2).
//
// Re-runs AllEquationsAboveCovers and compares what the pipeline actually produces, so passing IS
// reproduction -- the stronger claim than GuoYangEquations.m's stored-model comparison.
//
// ⚠ THIS BASE INCLUDES ITS `CRV` (paired) ENTRY, which was impossible before 2026-09-07. The
// helper now CONSTRUCTS the isomorphism for CRV pairs (tests/_crviso.m) instead of calling
// IsIsomorphic, which HANGS on them -- the genus-3 pair at 14_3 ran >50 min, the genus-5 one at
// 26_3 >1 h. The ambient weights are DERIVED (y's weight is half its own equation's degree), which
// tests/CRVStructure.m verifies against every stored entry.
//
// Keys with several stored entries are omitted: they cannot be matched unambiguously.
// The matrix in each cover_data value is a placeholder -- with manual_isomorphism false (the
// default) the helper never reads it.
//
// ✅ INVOLUTIONS CHECKED (2026-09-08). Guo-Yang publish, in THEIR coordinates:
//     w_2(x,y,z) = (-x,-y,-z)   w_3(x,y,z) = (x,y,-z)   w_7(x,y,z) = (x,-y,z)
// Our model is a different presentation, so those do NOT carry over. They were TRANSPORTED:
// psi := IsIsomorphic(our stored pair, Guo-Yang's pair), computed from the two EQUATIONS alone,
// and the matrix recorded here is psi^-1 . w_GY . psi.
// ⚠ WHY THIS IS NOT CIRCULAR: the involutions are Guo-Yang's (external) and psi comes from
// equations, never from the pipeline's own `ws`. The harness then checks that the PIPELINE's
// involution labelled w_m matches Guo-Yang's w_m, so a labelling error is detectable.
//
// ⚠ THREE TRAPS HERE, each of which cost a wrong conclusion or a rerun:
//   * construct_crv_isomorphism DECLINES on this base -- Guo-Yang's y has weight 3 (a genus-2
//     y-quotient) against our weight 2 (genus 1), so the two present the curve over DIFFERENT
//     intermediate quotients and there is no common base to take a Mobius map from. The general
//     IsIsomorphic is the fallback (208 s here).
//   * Inverse(psi) raises "Map has no inverse" on the map IsIsomorphic returns. IsInvertible
//     succeeds on the SAME map -- a representation issue, not a mathematical one.
//   * The composite psi^-1 . w_GY . psi comes back as ONE unreduced degree-39 representation, so
//     reading coefficients off it reports "not linear" although the MAP is linear. These matrices
//     were SOLVED FOR instead: on P(1,2,1,1) only y has weight 2, so a weight-respecting matrix
//     must send y -> c*y and act on (x,s,z) by a 3x3 block, and q1*L3-q3*L1, q1*L4-q4*L1 vanishing
//     on the curve are LINEAR conditions on that block's 9 coefficients. Kernel dimension came out
//     1 (so the block is unique up to scalar) and each result was certified by MAP EQUALITY
//     against the transported map, which is representation-independent.

function load_covers_and_ws_data_21_2()
    _<s> := PolynomialRing(Rationals());

    P3_1<x,y,s,z> := WeightedProjectiveSpace(Rationals(), [1,2,1,1]);

    cover_data := AssociativeArray();
    cover_data[{1,42}] := <HyperellipticCurve(Polynomial(Rationals(), [ -7/16, 0, -1/32, 0, 9/256 ])), DiagonalMatrix([1,1,1])>;   // genus 1
    cover_data[{1,21}] := <HyperellipticCurve(Polynomial(Rationals(), [ -7/256, 0, 31/128, 0, 9/256 ])), DiagonalMatrix([1,1,1])>;   // genus 1
    cover_data[{1,2,7,14}] := <HyperellipticCurve(Polynomial(Rationals(), [ 0, 81/4, -81/4 ])), DiagonalMatrix([1,1,1])>;   // genus 0
    cover_data[{1,6,7,42}] := <HyperellipticCurve(Polynomial(Rationals(), [ 0, -3 ])), DiagonalMatrix([1,1,1])>;   // genus 0
    cover_data[{1,3,7,21}] := <HyperellipticCurve(Polynomial(Rationals(), [ -3, 3 ])), DiagonalMatrix([1,1,1])>;   // genus 0
    cover_data[{1,2,3,6}] := <HyperellipticCurve(Polynomial(Rationals(), [ -4/9, 32/9, -71/9, 28/9 ])), DiagonalMatrix([1,1,1])>;   // genus 1
    cover_data[{1,2,21,42}] := <HyperellipticCurve(Polynomial(Rationals(), [ -7/16, 3/32, 81/256 ])), DiagonalMatrix([1,1,1])>;   // genus 0
    cover_data[{1,3,14,42}] := <HyperellipticCurve(Polynomial(Rationals(), [ 0, 7/4, -1/8, -9/64 ])), DiagonalMatrix([1,1,1])>;   // genus 1
    cover_data[{1,6,14,21}] := <HyperellipticCurve(Polynomial(Rationals(), [ 0, 7/64, 31/32, -9/64 ])), DiagonalMatrix([1,1,1])>;   // genus 1
    cover_data[{1}] := <Curve(P3_1, [ y^2 + 7/256*s^4 + 25/32*s^2*z^2 + 7/16*z^4, x^2 + 3*s^2 + 3*z^2 ]), DiagonalMatrix([1,1,1,1])>;   // genus 3, CRV pair

    ws_data := AssociativeArray();
    ws_data[{1}] := AssociativeArray();
    ws_data[{1}][2] := Matrix(4,4,[ 1,0,0,0,  0, 1,0,0,  0,0,-1,0,  0,0,0,-1 ]);
    ws_data[{1}][3] := Matrix(4,4,[ 1,0,0,0,  0,-1,0,0,  0,0,-1,0,  0,0,0, 1 ]);
    ws_data[{1}][7] := Matrix(4,4,[ 1,0,0,0,  0,-1,0,0,  0,0, 1,0,  0,0,0, 1 ]);
    return cover_data, ws_data;
end function;

procedure test_21_2()
    cover_data, ws_data := load_covers_and_ws_data_21_2();
    curves := GetHyperellipticCandidates();
    // ⚠ base_label := 7079 IS REQUIRED, and the reason is the V_4.
    // A CRV pair is built over a CHOSEN BASE, and the base decides which Klein four-group the pair
    // presents. The DEFAULT run builds W={1} over base 7078 (conic -x^2-3) and yields a pair whose
    // y-quotient is NOT isomorphic to the stored one and whose conic differs by a non-square -- a
    // different, equally valid V_4. Base 7079 (conic -3x^2-3) reproduces the committed entry
    // exactly. Same lesson as 26_3: the presentation is a choice, the curve is not.
    test_AllEquationsAboveCoversSingleCurve(21, 2, cover_data, ws_data, curves : base_label := 7079);
    return;
end procedure;

test_21_2();
