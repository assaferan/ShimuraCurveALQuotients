import "tests/BorcherdsProducts.m" : test_AllEquationsAboveCoversSingleCurve;

function load_covers_and_ws_data_14_1()
    _<x> := PolynomialRing(Rationals());

// EXTERNAL SOURCE: Gonzalez-Rotger, "Non-elliptic Shimura curves of genus one",
// J. Math. Soc. Japan 58 (2006) 927-948; Table 1 p.8 and the involution table directly below it.
// ⚠ THIS IS THE FIRST X0_ TEST BUILT ON GONZALEZ-ROTGER RATHER THAN GUO-YANG. The X0_ batches
// swept Guo-Yang's 43 bases, and GR's eleven genus-one bases are DISJOINT from those 43, so ten of
// the eleven had no re-derivation test at all -- the same blind spot as tests/OracleCoverage.m's,
// one level up: a sweep is blind to the bases it does not enumerate. (15_1 is the exception; it
// already had one.)
//
// ✅ INVOLUTIONS PUBLISHED, NOT DERIVED -- unusually comfortable for this repo. GR print the action
// of every Atkin-Lehner involution on their own model (p.8): w_{D.N}(x,y) = (x,-y) for all eleven
// curves, and for D = 14 the remaining generator is w_2(x,y) = (-x,y). So
//     w_2(x,y) = (-x, y)     w_7 = w_2*w_14 = (-x,-y)     w_14(x,y) = (x,-y)
// ⚠ NO TRANSPORT IS NEEDED HERE, and that is WHY 14_1 was chosen as the pilot: cover_data[{1}] is
// GR's equation VERBATIM, so their coordinates ARE the expected curve's coordinates and their
// involutions can be written down as matrices without carrying anything through an isomorphism.
// The bases with non-diagonal involutions (34_1, 10_7, and the four with N > 1) DO need the
// transport of [[gy-involution-transport]] and are deliberately not attempted here.
//
// ⚠ WHAT PINS THE LABELLING, since w_2 and w_14 both have genus-0 quotients and genus alone cannot
// separate them: the THIRD involution does. GR's model forces X/w_2 to be genus 0 (f is even in x,
// so u = x^2 descends), X/w_7 to be genus 1 (u = x^2, v = x*y gives v^2 = u*f, a cubic), and
// X/w_14 to be genus 0 (the x-line, P^1 by their Lemma 2.1). Our committed model has exactly that
// genus pattern 0/1/0, so a swap of w_2 and w_7 is already excluded before the pipeline runs -- and
// the quotient equations below make the check an equality, not just a genus count.
//
// ⚠ VERIFIED AGAINST THE COMMITTED MODEL BEFORE THIS TEST WAS WRITTEN, so that a red result here
// means the RE-DERIVATION disagrees, not that the expected curves were wrong:
//     W=[1]    proved Q-isomorphic to GR's quartic       (tests/GonzalezRotger.m PART 1)
//     W=[1,2]  conic class [7] both sides                (tests/GonzalezRotger.m PART 2)
//     W=[1,7]  IsIsomorphic TRUE against v^2 = u*f(u)    -- NOT covered by GonzalezRotger.m, whose
//              PART 2 only treats the even-quartic w_2-type quotient
//     W=[1,14] split, i.e. P^1 over Q                    (tests/GonzalezRotger.m PART 3)

    // verifying [Gonzalez-Rotger, Table 1, p. 8]
    // D = 14
    cover_data := AssociativeArray();
    // ⚠ The second component is only read when manual_isomorphism is set, which it is not here;
    // the helper decides by IsIsomorphic (genus > 0) or conic class (genus 0). It is the identity
    // rather than a made-up matrix so that nothing here looks like a verified coordinate change.
    id3 := IdentityMatrix(Rationals(), 3);
    cover_data[{1}]     := <HyperellipticCurve(-x^4 + 13*x^2 - 128),      id3>;
    cover_data[{1,2}]   := <HyperellipticCurve(-x^2 + 13*x - 128),        id3>;   // u = x^2
    cover_data[{1,7}]   := <HyperellipticCurve(x*(-x^2 + 13*x - 128)),    id3>;   // u = x^2, v = x*y
    cover_data[{1,14}]  := <HyperellipticCurve(x^2 - x),                  id3>;   // P^1: the split class

    ws_data := AssociativeArray();
    ws_data[{1}] := AssociativeArray();
    ws_data[{1}][2]  := DiagonalMatrix([-1, 1, 1]);
    ws_data[{1}][7]  := DiagonalMatrix([-1,-1, 1]);
    ws_data[{1}][14] := DiagonalMatrix([ 1,-1, 1]);
    // ⚠ No ws_data for the genus-0 keys: the helper's genus-0 branch decides by conic class and
    // exhibits no isomorphism to conjugate by, and it fails loudly rather than skip if asked.

    return cover_data, ws_data;
end function;

procedure test_14_1()
    cover_data, ws_data := load_covers_and_ws_data_14_1();
    curves := GetHyperellipticCandidates();

    test_AllEquationsAboveCoversSingleCurve(14, 1, cover_data, ws_data, curves);
    return;
end procedure;

test_14_1();
