import "tests/BorcherdsProducts.m" : test_AllEquationsAboveCoversSingleCurve;

function load_covers_and_ws_data_34_1()
    _<x> := PolynomialRing(Rationals());

// EXTERNAL SOURCE: Gonzalez-Rotger, "Non-elliptic Shimura curves of genus one",
// J. Math. Soc. Japan 58 (2006) 927-948; Table 1 p.8 and the involution table directly below it.
// Built on GONZALEZ-ROTGER, not Guo-Yang -- see tests/X0_14_1.m for why ten of GR's eleven
// genus-one bases had no re-derivation test until 2026-09-23.
//
// ✅ INVOLUTIONS PUBLISHED, NOT DERIVED (p.8):
//     w_34(x,y) = (x, -y)      w_17(x,y) = (-1/x, -y/x^2)      w_2 = w_17*w_34 = (-1/x, y/x^2)
// Lemma 3.2's I_0 cell for (34,1) is {w_34, w_17} with rows K_34 and K_17 -- subscripts matching,
// unlike the (21,1) cell, whose typo is documented in tests/_offline/X0_21_1.m.
//
// ⚠ 34_1 IS THE FIRST OF THESE WITH NON-DIAGONAL INVOLUTIONS, AND IT NEEDED NO TRANSPORT. That is
// worth stating because PLAN predicted otherwise. The transport problem of [[gy-involution-transport]]
// arises when the published involution must be carried into OUR model's coordinates. Here it must
// not: the helper conjugates by the isomorphism it finds and compares against ws_data in the
// EXPECTED curve's coordinates, and cover_data[{1}] is GR's equation VERBATIM. So GR's formulas are
// already in the right frame; the only work is writing them as matrices.
//
// ⚠ THE CONVENTION IS ROW-VECTOR TIMES MATRIX -- the helper forms
// `Vector(x)*ChangeRing(ws_ex[Q], Universe(x))` on the ambient coordinates (X:Y:Z) of the weighted
// projective space P(1,2,1), where the affine point is x = X/Z, y = Y/Z^2. So
//     (x,y) -> (-1/x, -y/x^2)   is   (X:Y:Z) -> (-Z : -Y : X)
//     (x,y) -> (-1/x,  y/x^2)   is   (X:Y:Z) -> (-Z :  Y : X)
// ✅ VERIFIED, NOT ASSERTED: all three matrices are automorphisms of GR's curve AND involutions,
// and they satisfy the labelled group law w_17 * w_34 = w_2, checked by IsIsomorphism and map
// equality before this file was written.
//
// ⚠ THE QUOTIENT EQUATIONS ARE DERIVED, so they were checked symbolically and then against the
// committed model. With u = x - 1/x (invariant under x -> -1/x):
//     w_17: v = y/x is invariant, and v^2 = f(x)/x^2 = -3u^2 + 26u - 59      [genus 0]
//     w_2 : s = x + 1/x is ANTI-invariant with s^2 = u^2 + 4, so t = (y/x)*s is invariant and
//           t^2 = (-3u^2 + 26u - 59)(u^2 + 4)                                [genus 1]
// Both identities verified in the function field; both quotients then verified against the
// committed model (conic class [2] for W=[1,17]; IsIsomorphic TRUE for W=[1,2]); W=[1,34] split.

    // verifying [Gonzalez-Rotger, Table 1, p. 8]
    // D = 34
    cover_data := AssociativeArray();
    // ⚠ The second component is read only when manual_isomorphism is set, which it is not here.
    id3 := IdentityMatrix(Rationals(), 3);
    cover_data[{1}]    := <HyperellipticCurve(-3*x^4 + 26*x^3 - 53*x^2 - 26*x - 3), id3>;
    cover_data[{1,17}] := <HyperellipticCurve(-3*x^2 + 26*x - 59),                  id3>;  // v = y/x
    cover_data[{1,2}]  := <HyperellipticCurve((-3*x^2 + 26*x - 59)*(x^2 + 4)),      id3>;  // t = (y/x)(x+1/x)
    cover_data[{1,34}] := <HyperellipticCurve(x^2 - x),                             id3>;  // P^1: split class

    ws_data := AssociativeArray();
    ws_data[{1}] := AssociativeArray();
    ws_data[{1}][34] := DiagonalMatrix(Rationals(), [1,-1,1]);
    ws_data[{1}][17] := Matrix(Rationals(), 3, 3, [0,0,1, 0,-1,0, -1,0,0]);
    ws_data[{1}][2]  := Matrix(Rationals(), 3, 3, [0,0,1, 0, 1,0, -1,0,0]);
    // ⚠ No ws_data for the genus-0 keys: the helper's genus-0 branch decides by conic class and
    // exhibits no isomorphism to conjugate by, and fails loudly rather than skip if asked.

    return cover_data, ws_data;
end function;

procedure test_34_1()
    cover_data, ws_data := load_covers_and_ws_data_34_1();
    curves := GetHyperellipticCandidates();

    test_AllEquationsAboveCoversSingleCurve(34, 1, cover_data, ws_data, curves);
    return;
end procedure;

test_34_1();
