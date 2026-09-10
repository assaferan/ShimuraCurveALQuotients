import "tests/BorcherdsProducts.m" : test_AllEquationsAboveCoversSingleCurve;

// tests/X0_14_3.m -- RE-DERIVATION test for X_0^14(3).
//
// Unlike tests/GuoYangEquations.m, which reads the STORED model, this re-runs
// AllEquationsAboveCovers and compares what the pipeline actually produces. Passing therefore
// IS reproduction, which is the stronger claim (see data/models/PROVENANCE.md).
//
// The expected curves below are the committed model's hyperelliptic cover entries, which are
// themselves validated against Guo-Yang's published equations. CRV (paired) entries and empty
// keys are omitted: the helper compares hyperelliptic covers, and a key with several stored
// entries cannot be matched unambiguously.
//
//
// ⚠ WHAT A PASS DOES AND DOES NOT PIN, measured 2026-09-07 by negative control on X0_14_3.
// The comparison is `IsIsomorphic`, so it pins each cover UP TO ISOMORPHISM -- which is the right
// notion for a model, but it is WEAKER on the genus-0 entries than it looks:
//   * perturbing the GENUS-2 entry at {1,3} (-184 -> -185) makes the test FAIL, as it should;
//   * perturbing a GENUS-0 conic's constant term does NOT -- two conics can be isomorphic with
//     different coefficients, so `IsIsomorphic` correctly still says yes.
// So the genus-0 cover entries are checked only for their isomorphism CLASS. If a conic's actual
// coefficients matter (they do for the CRV constructions, and for `15_4`'s twist), that has to be
// pinned elsewhere -- see tests/CRVFullCurve.m and tests/CRV_15_4.m.
// The second component of each cover_data value is unused here -- with manual_isomorphism false
// (the default) the helper calls IsIsomorphic, so the matrix is a placeholder. ws_data is left
// empty for the same reason: the helper skips involution checks for keys it does not find.
//
// ✅ INVOLUTIONS CHECKED (2026-09-08), and the W={1} CRV pair is compared too. Guo-Yang publish:
//     w_2(x,y,z) = (-x,y,z)   w_3(x,y,z) = (x,-y,-z)   w_14(x,y,z) = (x,-y,z)
// Those are in THEIR coordinates and were TRANSPORTED into ours:
// psi := construct_crv_isomorphism(our pair, theirs), computed from the two EQUATIONS alone, and
// the matrix recorded here is psi^-1 . w_GY . psi. NOT circular: the involutions are Guo-Yang's
// and psi comes from equations, never from the pipeline's own `ws`, which is what the harness
// then checks against them -- so a labelling error in the pipeline is detectable.
//
// ⚠ base_label := 5394 IS REQUIRED, and finding it is the whole point of this entry. A CRV pair is
// built over a chosen base cover, and THE BASE DECIDES WHICH V_4 THE PAIR PRESENTS. The default
// base (5383) gives a valid pair that is NOT the one Guo-Yang present, and the two are then not
// directly comparable -- construct_crv_isomorphism declines, and the general IsIsomorphic needs
// 6739 s (112 min) to confirm they are abstractly isomorphic while still yielding no usable
// coordinate change. Sweeping the candidate bases (7 s each, replaying only the pointless-conic
// step) finds 5394, whose pair IS Guo-Yang's V_4:
//     ours  y^2 = 4s^4 + 88s^2z^2 - 28z^4,   x^2 = -2s^2 - 9z^2
//     GY    y^2 = -7x^4 + 22x^2 + 1,         z^2 = -9x^2 - 2
// i.e. theirs scaled by 4 with the two base coordinates exchanged. The constructor then finds the
// isomorphism immediately. ⇒ WHEN A PAIR WILL NOT MATCH, SWEEP THE BASE BEFORE CONCLUDING
// ANYTHING ABOUT THE CURVE -- same lesson as 26_3's base_label := 8103.
//
// ⚠ The matrices were SOLVED FOR, not read off: the composite is an unreduced representation, so
// inspecting its coefficients reports "not linear" although the MAP is linear. On P(1,2,1,1) only
// y has weight 2, so a weight-respecting matrix sends y -> c*y and acts on (x,s,z) by a 3x3 block;
// q1*L3-q3*L1 and q1*L4-q4*L1 vanishing on the curve are LINEAR in that block's 9 coefficients.
// Kernel dimension 1, and each result certified by MAP EQUALITY. Script: tests/_gyinvol_crv.m.

function load_covers_and_ws_data_14_3()
    _<s> := PolynomialRing(Rationals());

    P3_1<x,y,s,z> := WeightedProjectiveSpace(Rationals(), [1,2,1,1]);

    cover_data := AssociativeArray();
    cover_data[{1}] := <Curve(P3_1, [ y^2 - 4*s^4 - 88*s^2*z^2 + 28*z^4, x^2 + 2*s^2 + 9*z^2 ]), DiagonalMatrix([1,1,1,1])>;   // genus 3, CRV pair -- the base_label 5394 presentation
    cover_data[{1,2,7,14}] := <HyperellipticCurve(Polynomial(Rationals(), [ -11, 104, -128 ])), DiagonalMatrix([1,1,1])>;   // genus 0
    cover_data[{1,42}] := <HyperellipticCurve(Polynomial(Rationals(), [ -28, 0, 88, 0, 4 ])), DiagonalMatrix([1,1,1])>;   // genus 1
    cover_data[{1,3,7,21}] := <HyperellipticCurve(Polynomial(Rationals(), [ -176, 2368, -3776, 1024 ])), DiagonalMatrix([1,1,1])>;   // genus 1
    cover_data[{1,3}] := <HyperellipticCurve(Polynomial(Rationals(), [ 63, 0, -184, 0, -53, 0, -2 ])), DiagonalMatrix([1,1,1])>;   // genus 2
    cover_data[{1,7}] := <HyperellipticCurve(Polynomial(Rationals(), [ 768, 2304, 4160, 45440/9, 4160, 2304, 768 ])), DiagonalMatrix([1,1,1])>;   // genus 2
    cover_data[{1,6}] := <HyperellipticCurve(Polynomial(Rationals(), [ 3087/2, 0, 577/2, 0, 17/2, 0, -1/2 ])), DiagonalMatrix([1,1,1])>;   // genus 2
    cover_data[{1,2,3,6}] := <HyperellipticCurve(Polynomial(Rationals(), [ -11, 236, -1420, 1952, -512 ])), DiagonalMatrix([1,1,1])>;   // genus 1
    cover_data[{1,2,21,42}] := <HyperellipticCurve(Polynomial(Rationals(), [ 64, -768, 256 ])), DiagonalMatrix([1,1,1])>;   // genus 0
    cover_data[{1,6,14,21}] := <HyperellipticCurve(Polynomial(Rationals(), [ -11, 16 ])), DiagonalMatrix([1,1,1])>;   // genus 0
    cover_data[{1,3,14,42}] := <HyperellipticCurve(Polynomial(Rationals(), [ 1, -8 ])), DiagonalMatrix([1,1,1])>;   // genus 0
    cover_data[{1,21}] := <HyperellipticCurve(Polynomial(Rationals(), [ -343, 0, -26, 0, 1 ])), DiagonalMatrix([1,1,1])>;   // genus 1
    cover_data[{1,6,7,42}] := <HyperellipticCurve(Polynomial(Rationals(), [ 64, -1280, 6400, -2048 ])), DiagonalMatrix([1,1,1])>;   // genus 1

    ws_data := AssociativeArray();
    ws_data[{1}] := AssociativeArray();
    ws_data[{1}][2]  := Matrix(4,4,[ 1,0,0,0,  0, 1,0,0,  0,0, 1,0,  0,0,0,-1 ]);
    ws_data[{1}][3]  := Matrix(4,4,[ 1,0,0,0,  0,-1,0,0,  0,0,-1,0,  0,0,0,-1 ]);
    ws_data[{1}][14] := Matrix(4,4,[ 1,0,0,0,  0,-1,0,0,  0,0, 1,0,  0,0,0, 1 ]);
    return cover_data, ws_data;
end function;

procedure test_14_3()
    cover_data, ws_data := load_covers_and_ws_data_14_3();
    curves := GetHyperellipticCandidates();
    test_AllEquationsAboveCoversSingleCurve(14, 3, cover_data, ws_data, curves : base_label := 5394);
    return;
end procedure;

test_14_3();
