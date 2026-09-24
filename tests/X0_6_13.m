import "tests/BorcherdsProducts.m" : test_AllEquationsAboveCoversSingleCurve;

function load_covers_and_ws_data_6_13()
    _<x> := PolynomialRing(Rationals());

// EXTERNAL SOURCE: Gonzalez-Rotger, "Non-elliptic Shimura curves of genus one",
// J. Math. Soc. Japan 58 (2006) 927-948; Table 1 p.8 and the involution table below it.
// Built on GONZALEZ-ROTGER, not Guo-Yang -- see tests/X0_14_1.m for the coverage gap this closes.
//
// ⚠⚠ THE SECOND OF THE TWO BASES NEEDING A LIST IN cover_data. W={1} comes out over TWO bases
// (1532 and 1533) whose quartics are NOT isomorphic as CrvHyp -- different degree-2 divisor classes
// of the same curve. Base 1532 gives GR's published curve.
// ⇒ ENTRY 1 IS THE ORACLE (GR's quartic), ENTRY 2 A DRIFT ALTERNATIVE (our committed entry 2, which
// nothing external vouches for). The helper requires every produced cover to match one of them and
// at least one to match ENTRY 1.
// ⚠ Note the asymmetry with tests/X0_6_5.m: there the PUBLISHED curve is our committed entry 2 and
// the alternative is entry 1; here it is the other way round. tests/GonzalezRotger.m records the
// pair as NO_EXHIBITED_ISO = {<6,5,1>, <6,13,2>} for exactly that reason. Do not copy one file's
// ordering into the other.
// ⚠ PINNING base_label DOES NOT WORK, measured: it is threaded only into
// EquationsAbovePointlessConics and EquationsByRebase, while W={1} comes from the earlier step.
//
// ✅ INVOLUTIONS PUBLISHED (p.8): w_78(x,y) = (x,-y), w_2(x,y) = (-x,y),
// w_3(x,y) = (64/x, 64y/x^2). Lemma 3.2 gives I_0(X_0(6,13)) = {w_78, w_2, w_3, w_13}, so those
// four have genus-0 quotients and w_6, w_26, w_39 have genus-1 quotients.
// ⚠ COMPOSE CAREFULLY: w_3(w_2(x,y)) = w_3(-x,y) = (-64/x, 64y/x^2) = w_6.
//
// ⚠ QUOTIENT EQUATIONS DERIVED, identities checked in the function field:
//     u = x + 64/x : f/x^2 = -u^2 + 13       (w_3  : v = y/x)
//     w = x - 64/x : f/x^2 = -w^2 - 243      (w_13 : v = y/x)
// Both, plus W=[1,2] (u = x^2), verified against the committed model: classes [3], [], [3] on both
// sides.
// ⚠ The genus-1 quotient keys (w_6, w_26, w_39) are deliberately absent for the same reason as in
// tests/X0_6_7.m and tests/X0_10_3.m -- same Jacobian, inequivalent quartic model; the helper's
// second pass drift-checks them regardless.

    // verifying [Gonzalez-Rotger, Table 1, p. 8]
    // D = 6, N = 13
    cover_data := AssociativeArray();
    id3 := IdentityMatrix(Rationals(), 3);   // read only under manual_isomorphism, which is off
    cover_data[{1}] := <[* HyperellipticCurve(-x^4 - 115*x^2 - 4096),          // 1: GR, base 1532
                           HyperellipticCurve(-10647/7396*x^4 + 507/3698*x^3
                                              - 20449/7396*x^2 - 507/3698*x
                                              - 10647/7396) *],               // 2: ours, base 1533
                        id3>;
    cover_data[{1,2}]  := <HyperellipticCurve(-x^2 - 115*x - 4096), id3>;  // u = x^2
    cover_data[{1,3}]  := <HyperellipticCurve(-x^2 + 13),           id3>;  // u = x+64/x, v = y/x
    cover_data[{1,13}] := <HyperellipticCurve(-x^2 - 243),          id3>;  // w = x-64/x, v = y/x
    cover_data[{1,78}] := <HyperellipticCurve(x^2 - x),             id3>;  // P^1

    M78 := DiagonalMatrix(Rationals(), [1,-1,1]);                  // (x, -y)
    M2  := DiagonalMatrix(Rationals(), [-1,1,1]);                  // (-x, y)
    M3  := Matrix(Rationals(), 3, 3, [0,0,1, 0,64,0, 64,0,0]);     // (64/x, 64y/x^2)

    ws_data := AssociativeArray();
    ws_data[{1}] := AssociativeArray();
    ws_data[{1}][78] := M78;
    ws_data[{1}][2]  := M2;
    ws_data[{1}][3]  := M3;
    ws_data[{1}][39] := M2 * M78;         // w_2 * w_78 = w_{156/4}  = w_39
    ws_data[{1}][26] := M3 * M78;         // w_3 * w_78 = w_{234/9}  = w_26
    ws_data[{1}][6]  := M2 * M3;          // w_2 * w_3  = w_6
    ws_data[{1}][13] := M2 * M3 * M78;    // w_6 * w_78 = w_{468/36} = w_13
    // ⚠ No ws_data for the genus-0 keys: the helper's genus-0 branch exhibits no isomorphism to
    // conjugate by, and fails loudly rather than skip if asked.

    return cover_data, ws_data;
end function;

procedure test_6_13()
    cover_data, ws_data := load_covers_and_ws_data_6_13();
    curves := GetHyperellipticCandidates();

    test_AllEquationsAboveCoversSingleCurve(6, 13, cover_data, ws_data, curves);
    return;
end procedure;

test_6_13();
