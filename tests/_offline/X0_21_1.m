import "tests/BorcherdsProducts.m" : test_AllEquationsAboveCoversSingleCurve;

function load_covers_and_ws_data_21_1()
    _<x> := PolynomialRing(Rationals());

// EXTERNAL SOURCE: Gonzalez-Rotger, "Non-elliptic Shimura curves of genus one",
// J. Math. Soc. Japan 58 (2006) 927-948; Table 1 p.8 and the involution table directly below it.
// Built on GONZALEZ-ROTGER, not Guo-Yang -- see tests/X0_14_1.m for why ten of GR's eleven
// genus-one bases had no re-derivation test until 2026-09-23.
//
// ⚠⚠ THIS TEST IS IN _offline/ FOR ONE REASON: IT REQUIRES HMFIT=1.
//
//     HMFIT=1 magma -b filename:=tests/_offline/X0_21_1.m run_tests.m < /dev/null     # 44 s
//
// Under DEFAULT flags the pipeline cannot re-derive 21_1 at all. It dies before any comparison,
// in find_signs_hauptmodul:
//     "no choice of signs satisfies s/scale + stilde/scale_tilde = 1 at discriminant(s)
//      [-15, -43, -51, -67]"
// That is not a defect in the model and not a defect in this test -- it is the documented reason
// 21_1 was built under HMFIT in the first place (HMFIT solves for (scale, scale_tilde) against
// every rational CM point instead of reading them off the normalising points). HMFIT is env-gated
// and OFF by default, so a CI-visible copy of this test would be permanently red.
// ⇒ **21_1's committed model is CORRECT BUT NOT REPRODUCIBLE UNDER DEFAULT FLAGS.** That is a
// property of the pipeline's default normalisation, not of the model: with HMFIT=1 this test
// re-derives all four covers and matches every published involution.
// ⚠ Do NOT "fix" this by deleting the flag requirement or by making the test skip itself when the
// flag is absent -- a test that silently skips is the vacuity failure mode this repo has been
// bitten by. If HMFIT ever becomes the default, move this file into tests/.
//
// ✅ INVOLUTIONS PUBLISHED, NOT DERIVED. GR print the action on their own model (p.8):
// w_{D.N}(x,y) = (x,-y) for all eleven curves, and for D = 21 the remaining generator is
// w_7(x,y) = (-x,y).  Hence
//     w_7(x,y) = (-x, y)     w_3 = w_7*w_21 = (-x,-y)     w_21(x,y) = (x,-y)
//
// ⚠⚠ A PUBLISHED TYPO, AND IT IS IN THE INVOLUTION SET FOR THIS VERY BASE -- read this before
// "correcting" w_7 to w_3. Lemma 3.2 (p.7) prints the I_0 set for (21,1) as {w_21, w_3}. But the
// two rows of that same cell are K_21 and K_7, NOT K_3, and every other cell in the lemma has
// matching subscripts. The set is {w_21, w_7}. Three independent things agree:
//   * the paper's own K_7 row, inside the inconsistent cell;
//   * the p.8 involution table, which gives w_7(x,y) = (-x,y) for D = 21;
//   * our committed model: genus 0 at W=[1,7] and genus 1 at W=[1,3] -- and a genus-0 quotient is
//     exactly what membership of I_0 means.
// GR's own p.8 Jacobian table corroborates: it records Jac(X_0(21,1)/<u.w>) = 21A6, which is genus
// one, and u.w is w_3 precisely when u = w_7.
//
// ⚠ NO TRANSPORT IS NEEDED: cover_data[{1}] is GR's equation VERBATIM, so their coordinates ARE the
// expected curve's coordinates and the involutions are diagonal matrices. The bases whose
// involutions are non-diagonal (34_1, 10_7, and the four at N > 1) need the transport of
// [[gy-involution-transport]] and are NOT in this group.
//
// ⚠ WHAT PINS THE LABELLING, since w_7 and w_21 both have genus-0 quotients and genus alone cannot
// separate them: the THIRD involution does. GR's quartic is EVEN in x, so u = x^2 descends and
// X/w_7 is the conic y^2 = f(u); X/w_3 takes u = x^2 and v = x*y, giving v^2 = u*f(u), which is
// genus 1; and X/w_21 is the x-line, P^1 over Q by their Lemma 2.1.
//
// ⚠ VERIFIED AGAINST THE COMMITTED MODEL BEFORE THIS TEST WAS WRITTEN, so a red result here means
// the RE-DERIVATION disagrees, not that the expected curves are wrong:
//     W=[1]     proved Q-isomorphic to GR's quartic       (tests/GonzalezRotger.m PART 1)
//     W=[1,7]   conic class [3] both sides                (tests/GonzalezRotger.m PART 2)
//     W=[1,3]   IsIsomorphic TRUE against v^2 = u*f(u)    -- NOT covered by GonzalezRotger.m,
//               whose PART 2 only treats the even-quartic w_7-type quotient
//     W=[1,21]  split, i.e. P^1 over Q                    (tests/GonzalezRotger.m PART 3)

    // verifying [Gonzalez-Rotger, Table 1, p. 8]
    // D = 21
    cover_data := AssociativeArray();
    // ⚠ The second component is read only when manual_isomorphism is set, which it is not here;
    // the helper decides by IsIsomorphic (genus > 0) or conic class (genus 0). Identity rather than
    // a made-up matrix, so nothing here looks like a verified coordinate change.
    id3 := IdentityMatrix(Rationals(), 3);
    cover_data[{1}]    := <HyperellipticCurve(-7*x^4 + 94*x^2 - 343),    id3>;
    cover_data[{1,7}]  := <HyperellipticCurve(-7*x^2 + 94*x - 343),      id3>;  // u = x^2
    cover_data[{1,3}]  := <HyperellipticCurve(x*(-7*x^2 + 94*x - 343)),  id3>;  // u = x^2, v = x*y
    cover_data[{1,21}] := <HyperellipticCurve(x^2 - x),                  id3>;  // P^1: split class

    ws_data := AssociativeArray();
    ws_data[{1}] := AssociativeArray();
    ws_data[{1}][7]  := DiagonalMatrix([-1, 1, 1]);
    ws_data[{1}][3]  := DiagonalMatrix([-1,-1, 1]);
    ws_data[{1}][21] := DiagonalMatrix([ 1,-1, 1]);
    // ⚠ No ws_data for the genus-0 keys: the helper's genus-0 branch decides by conic class and
    // exhibits no isomorphism to conjugate by, and fails loudly rather than skip if asked.

    return cover_data, ws_data;
end function;

procedure test_21_1()
    cover_data, ws_data := load_covers_and_ws_data_21_1();
    curves := GetHyperellipticCandidates();

    test_AllEquationsAboveCoversSingleCurve(21, 1, cover_data, ws_data, curves);
    return;
end procedure;

test_21_1();
