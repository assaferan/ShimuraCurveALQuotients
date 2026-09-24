import "tests/BorcherdsProducts.m" : test_AllEquationsAboveCoversSingleCurve;

function load_covers_and_ws_data_10_3()
    _<x> := PolynomialRing(Rationals());

// EXTERNAL SOURCE: Gonzalez-Rotger, "Non-elliptic Shimura curves of genus one",
// J. Math. Soc. Japan 58 (2006) 927-948; Table 1 p.8 and the involution table below it.
// Built on GONZALEZ-ROTGER, not Guo-Yang -- see tests/X0_14_1.m for the coverage gap this closes.
//
// ⚠ 10_3 IS THE BASE WHOSE W=[1,2] DRIFT THE GONZALEZ-ROTGER ORACLE CAUGHT after three internal
// consistency tests had passed it for a week -- all three committed entries were consistently
// wrong, and only an external oracle could arbitrate. This test adds the re-derivation half:
// the oracle says the committed curves are right, this says the pipeline still produces them.
//
// ✅ INVOLUTIONS PUBLISHED (p.8): w_30(x,y) = (x,-y), w_2(x,y) = (-x,y),
// w_3(x,y) = (4/x, 4y/x^2). Lemma 3.2 gives I_0(X_0(10,3)) = {w_30, w_2, w_3, w_5}, so those four
// have genus-0 quotients and w_6, w_10, w_15 have genus-1 quotients -- as our committed model shows.
//
// ⚠ COMPOSE CAREFULLY. w_3(w_2(x,y)) = w_3(-x,y) = (-4/x, 4y/x^2) = w_6. Getting this wrong at the
// sister base 6_7 produced a quotient Jacobian that disagreed with the model, which reads as
// "the model is wrong" when only the LABEL was wrong. Correctly composed, every key matches a
// committed entry:
//     w_2 class [5]    w_3 class [2,5]   w_5 class [2]                 (genus 0)
//     w_15 Jac 30a1    w_10 Jac 30a4     w_6 Jac 30a5                  (genus 1)
// The four Mobius maps x, -x, 4/x, -4/x, each lifting by y -> +-y, give all eight elements of
// W_{10,3} = {1,2,3,5,6,10,15,30}.
//
// ⚠ NO TRANSPORT IS NEEDED despite the non-diagonal maps: cover_data[{1}] is GR's equation
// VERBATIM. Convention is row-vector times matrix on (X:Y:Z) of P(1,2,1), x = X/Z, y = Y/Z^2.
// ✅ All seven non-trivial matrices verified as automorphisms AND involutions before this was
// written; the derived four come from w_m*w_n = w_{mn/gcd(m,n)^2}, so their labelling is forced.
//
// ⚠ QUOTIENT EQUATIONS ARE DERIVED; the identities were checked in the function field:
//     u = x + 4/x : f/x^2 = -2u^2 + 5      and   (x - 4/x)^2 = u^2 - 16
//     w = x - 4/x : f/x^2 = -2w^2 - 27     and   (x + 4/x)^2 = w^2 + 16
//
// ⚠ THE GENUS-1 KEYS EACH HOLD SEVERAL COMMITTED ENTRIES and only one of each is isomorphic to the
// curve derived from GR's model -- the others are inequivalent quartic models of the SAME curve
// (identical Jacobians), i.e. different degree-2 divisor classes. That is the torsor phenomenon
// documented in tests/GonzalezRotger.m, not a defect. It matters here only if the PIPELINE emits
// more than one cover at such a key; see the note in PLAN item 6 on the helper's one-curve-per-key
// contract.

    // verifying [Gonzalez-Rotger, Table 1, p. 8]
    // D = 10, N = 3
    cover_data := AssociativeArray();
    id3 := IdentityMatrix(Rationals(), 3);   // read only under manual_isomorphism, which is off
    cover_data[{1}]    := <HyperellipticCurve(-2*x^4 - 11*x^2 - 32),     id3>;
    cover_data[{1,2}]  := <HyperellipticCurve(-2*x^2 - 11*x - 32),       id3>;  // u=x^2,   v=y
    cover_data[{1,3}]  := <HyperellipticCurve(-2*x^2 + 5),               id3>;  // u=x+4/x, v=y/x
    cover_data[{1,5}]  := <HyperellipticCurve(-2*x^2 - 27),              id3>;  // w=x-4/x, v=y/x
    // ⚠ W=[1,6], W=[1,10] and W=[1,15] -- the three GENUS-1 quotient keys -- are DELIBERATELY
    // ABSENT, for the same reason W=[1,14] is absent from tests/X0_6_7.m. Each comes out over a
    // single base whose curve has the SAME Jacobian as the quotient derived from GR's model
    // (30a5, 30a4, 30a1 on both sides) but is not isomorphic to it as a CrvHyp: an inequivalent
    // quartic model of the same curve, i.e. a different degree-2 divisor class. That is the torsor
    // phenomenon documented in tests/GonzalezRotger.m, NOT a defect, and listing these keys would
    // make the test red for a reason that is not one.
    // ⚠ Nothing is lost by omitting them: the helper's SECOND pass compares every committed cover
    // key against the produced covers regardless, so they are still drift-checked. What is not
    // available for them is an ORACLE claim, because GR's published data pins the quotient's
    // Jacobian but not which degree-2 model of it the pipeline should choose.
    cover_data[{1,30}] := <HyperellipticCurve(x^2 - x),                  id3>;  // P^1

    M30 := DiagonalMatrix(Rationals(), [1,-1,1]);               // (x, -y)
    M2  := DiagonalMatrix(Rationals(), [-1,1,1]);               // (-x, y)
    M3  := Matrix(Rationals(), 3, 3, [0,0,1, 0,4,0, 4,0,0]);    // (4/x, 4y/x^2)

    ws_data := AssociativeArray();
    ws_data[{1}] := AssociativeArray();
    ws_data[{1}][30] := M30;
    ws_data[{1}][2]  := M2;
    ws_data[{1}][3]  := M3;
    ws_data[{1}][15] := M2 * M30;        // w_2 * w_30 = w_{60/4}  = w_15
    ws_data[{1}][10] := M3 * M30;        // w_3 * w_30 = w_{90/9}  = w_10
    ws_data[{1}][6]  := M2 * M3;         // w_2 * w_3  = w_6
    ws_data[{1}][5]  := M2 * M3 * M30;   // w_6 * w_30 = w_{180/36} = w_5
    // ⚠ No ws_data for the genus-0 keys: the helper's genus-0 branch exhibits no isomorphism to
    // conjugate by, and fails loudly rather than skip if asked.

    return cover_data, ws_data;
end function;

procedure test_10_3()
    cover_data, ws_data := load_covers_and_ws_data_10_3();
    curves := GetHyperellipticCandidates();

    test_AllEquationsAboveCoversSingleCurve(10, 3, cover_data, ws_data, curves);
    return;
end procedure;

test_10_3();
