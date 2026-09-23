import "tests/BorcherdsProducts.m" : test_AllEquationsAboveCoversSingleCurve;

function load_covers_and_ws_data_6_5()
    _<x> := PolynomialRing(Rationals());

// EXTERNAL SOURCE: Gonzalez-Rotger, "Non-elliptic Shimura curves of genus one",
// J. Math. Soc. Japan 58 (2006) 927-948; Table 1 p.8 and the involution table below it.
// Built on GONZALEZ-ROTGER, not Guo-Yang -- see tests/X0_14_1.m for the coverage gap this closes.
//
// ⚠⚠ THIS IS ONE OF THE TWO BASES THAT NEEDED cover_data TO LIST SEVERAL ACCEPTABLE CURVES.
// W={1} comes out over TWO bases (1483 and 1484) and the two quartics are NOT isomorphic as
// CrvHyp: they are different degree-2 divisor classes of the SAME curve -- a torsor under E(Q),
// which is non-trivial here. Base 1484 gives GR's published curve; base 1483 gives the other.
// ⇒ ENTRY 1 IS THE ORACLE (GR's published quartic) and ENTRY 2 IS A DRIFT ALTERNATIVE: our own
// committed entry 1, for which nothing external vouches. The helper requires every produced cover
// to match ONE of them and at least one to match ENTRY 1, so the published curve cannot quietly
// stop being consulted.
// ⚠ PINNING base_label DOES NOT WORK HERE and was measured, not assumed: base_label is threaded
// only into EquationsAbovePointlessConics and EquationsByRebase, while the W={1} cover comes from
// the earlier main step, so a run pinned at 1484 still emits BOTH bases.
// ⚠ tests/GonzalezRotger.m records this pair as NO_EXHIBITED_ISO = {<6,5,1>, <6,13,2>}. It is NOT
// a defect list: the two entries have identical everywhere-local solubility profiles, consistent
// with two models of one curve.
//
// ✅ INVOLUTIONS PUBLISHED (p.8): w_30(x,y) = (x,-y), w_2(x,y) = (-x,y),
// w_6(x,y) = (32/x, 32y/x^2). Lemma 3.2 gives I_0(X_0(6,5)) = {w_30, w_2, w_6, w_10}, so those four
// have genus-0 quotients and w_3, w_5, w_15 have genus-1 quotients.
// ⚠ COMPOSE CAREFULLY -- at the sister base 6_7 a wrong composition produced a quotient Jacobian
// that disagreed with the model and read as "the model is wrong" when only the label was.
// w_6(w_2(x,y)) = w_6(-x,y) = (-32/x, 32y/x^2) = w_3.
//
// ⚠ QUOTIENT EQUATIONS DERIVED, identities checked in the function field:
//     u = x + 32/x : f/x^2 = -u^2 + 125      (w_6  : v = y/x)
//     w = x - 32/x : f/x^2 = -w^2 - 3        (w_10 : v = y/x)
// Both genus-0 quotients, plus W=[1,2] (u = x^2), verified against the committed model: classes
// [3], [], [3] on both sides.
// ⚠ The three GENUS-1 quotient keys (w_3, w_5, w_15) are deliberately absent for the same reason as
// in tests/X0_6_7.m and tests/X0_10_3.m -- same Jacobian, inequivalent quartic model. The helper's
// second pass drift-checks them anyway; only an ORACLE claim is unavailable.

    // verifying [Gonzalez-Rotger, Table 1, p. 8]
    // D = 6, N = 5
    cover_data := AssociativeArray();
    id3 := IdentityMatrix(Rationals(), 3);   // read only under manual_isomorphism, which is off
    cover_data[{1}] := <[* HyperellipticCurve(-x^4 + 61*x^2 - 1024),          // 1: GR, base 1484
                           HyperellipticCurve(-64375/16384*x^4 + 15625/2048*x^3
                                              + 29375/8192*x^2 - 15625/2048*x
                                              - 64375/16384) *],             // 2: ours, base 1483
                        id3>;
    cover_data[{1,2}]  := <HyperellipticCurve(-x^2 + 61*x - 1024), id3>;   // u = x^2
    cover_data[{1,6}]  := <HyperellipticCurve(-x^2 + 125),         id3>;   // u = x+32/x, v = y/x
    cover_data[{1,10}] := <HyperellipticCurve(-x^2 - 3),           id3>;   // w = x-32/x, v = y/x
    cover_data[{1,30}] := <HyperellipticCurve(x^2 - x),            id3>;   // P^1

    M30 := DiagonalMatrix(Rationals(), [1,-1,1]);                  // (x, -y)
    M2  := DiagonalMatrix(Rationals(), [-1,1,1]);                  // (-x, y)
    M6  := Matrix(Rationals(), 3, 3, [0,0,1, 0,32,0, 32,0,0]);     // (32/x, 32y/x^2)

    ws_data := AssociativeArray();
    ws_data[{1}] := AssociativeArray();
    ws_data[{1}][30] := M30;
    ws_data[{1}][2]  := M2;
    ws_data[{1}][6]  := M6;
    ws_data[{1}][15] := M2 * M30;         // w_2 * w_30 = w_{60/4}   = w_15
    ws_data[{1}][5]  := M6 * M30;         // w_6 * w_30 = w_{180/36} = w_5
    ws_data[{1}][3]  := M2 * M6;          // w_2 * w_6  = w_{12/4}   = w_3
    ws_data[{1}][10] := M2 * M6 * M30;    // w_3 * w_30 = w_{90/9}   = w_10
    // ⚠ No ws_data for the genus-0 keys: the helper's genus-0 branch exhibits no isomorphism to
    // conjugate by, and fails loudly rather than skip if asked.

    return cover_data, ws_data;
end function;

procedure test_6_5()
    cover_data, ws_data := load_covers_and_ws_data_6_5();
    curves := GetHyperellipticCandidates();

    test_AllEquationsAboveCoversSingleCurve(6, 5, cover_data, ws_data, curves);
    return;
end procedure;

test_6_5();
