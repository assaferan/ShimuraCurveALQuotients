import "tests/BorcherdsProducts.m" : test_AllEquationsAboveCoversSingleCurve;

function load_covers_and_ws_data_46_1()
    _<x> := PolynomialRing(Rationals());

// EXTERNAL SOURCE: Gonzalez-Rotger, "Non-elliptic Shimura curves of genus one",
// J. Math. Soc. Japan 58 (2006) 927-948; Table 1 p.8 and the involution table directly below it.
// Built on GONZALEZ-ROTGER, not Guo-Yang -- see tests/X0_14_1.m for why ten of GR's eleven
// genus-one bases had no re-derivation test until 2026-09-23.
//
// ✅ INVOLUTIONS PUBLISHED, NOT DERIVED. GR print the action on their own model (p.8):
// w_{D.N}(x,y) = (x,-y) for all eleven curves, and for D = 46 the remaining generator is
// w_2(x,y) = (-x,y).  Hence
//     w_2(x,y) = (-x, y)     w_23 = w_2*w_46 = (-x,-y)     w_46(x,y) = (x,-y)
// Lemma 3.2's I_0 cell for (46,1) is {w_46, w_2}, and its rows are K_46 and K_2 -- subscripts
// matching, unlike the (21,1) cell, whose typo is documented in tests/X0_21_1.m.
//
// ⚠ NO TRANSPORT IS NEEDED: cover_data[{1}] is GR's equation VERBATIM, so their coordinates ARE the
// expected curve's coordinates and the involutions are diagonal matrices.
//
// ⚠ WHAT PINS THE LABELLING, since w_2 and w_46 both have genus-0 quotients and genus alone cannot
// separate them: the THIRD involution does. GR's quartic is EVEN in x, so u = x^2 descends and
// X/w_2 is the conic y^2 = f(u); X/w_23 takes u = x^2 and v = x*y, giving v^2 = u*f(u), which is
// genus 1; and X/w_46 is the x-line, P^1 over Q by their Lemma 2.1.
//
// ⚠ VERIFIED AGAINST THE COMMITTED MODEL BEFORE THIS TEST WAS WRITTEN, so a red result here means
// the RE-DERIVATION disagrees, not that the expected curves are wrong:
//     W=[1]     proved Q-isomorphic to GR's quartic       (tests/GonzalezRotger.m PART 1)
//     W=[1,2]   conic class [23] both sides               (tests/GonzalezRotger.m PART 2)
//     W=[1,23]  IsIsomorphic TRUE against v^2 = u*f(u)    -- NOT covered by GonzalezRotger.m,
//               whose PART 2 only treats the even-quartic w_2-type quotient
//     W=[1,46]  split, i.e. P^1 over Q                    (tests/GonzalezRotger.m PART 3)

    // verifying [Gonzalez-Rotger, Table 1, p. 8]
    // D = 46
    cover_data := AssociativeArray();
    // ⚠ The second component is read only when manual_isomorphism is set, which it is not here;
    // the helper decides by IsIsomorphic (genus > 0) or conic class (genus 0). Identity rather than
    // a made-up matrix, so nothing here looks like a verified coordinate change.
    id3 := IdentityMatrix(Rationals(), 3);
    cover_data[{1}]    := <HyperellipticCurve(-x^4 + 45*x^2 - 512),    id3>;
    cover_data[{1,2}]  := <HyperellipticCurve(-x^2 + 45*x - 512),      id3>;  // u = x^2
    cover_data[{1,23}] := <HyperellipticCurve(x*(-x^2 + 45*x - 512)),  id3>;  // u = x^2, v = x*y
    cover_data[{1,46}] := <HyperellipticCurve(x^2 - x),                id3>;  // P^1: split class

    ws_data := AssociativeArray();
    ws_data[{1}] := AssociativeArray();
    ws_data[{1}][2]  := DiagonalMatrix([-1, 1, 1]);
    ws_data[{1}][23] := DiagonalMatrix([-1,-1, 1]);
    ws_data[{1}][46] := DiagonalMatrix([ 1,-1, 1]);
    // ⚠ No ws_data for the genus-0 keys: the helper's genus-0 branch decides by conic class and
    // exhibits no isomorphism to conjugate by, and fails loudly rather than skip if asked.

    return cover_data, ws_data;
end function;

procedure test_46_1()
    cover_data, ws_data := load_covers_and_ws_data_46_1();
    curves := GetHyperellipticCandidates();

    test_AllEquationsAboveCoversSingleCurve(46, 1, cover_data, ws_data, curves);
    return;
end procedure;

test_46_1();
