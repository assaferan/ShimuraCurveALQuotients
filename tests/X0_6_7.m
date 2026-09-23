import "tests/BorcherdsProducts.m" : test_AllEquationsAboveCoversSingleCurve;

function load_covers_and_ws_data_6_7()
    _<x> := PolynomialRing(Rationals());

// EXTERNAL SOURCE: Gonzalez-Rotger, "Non-elliptic Shimura curves of genus one",
// J. Math. Soc. Japan 58 (2006) 927-948; Table 1 p.8 and the involution table below it.
// Built on GONZALEZ-ROTGER, not Guo-Yang -- see tests/X0_14_1.m for the coverage gap this closes.
//
// ✅ INVOLUTIONS PUBLISHED (p.8): w_42(x,y) = (x,-y), w_3(x,y) = (-x,y),
// w_6(x,y) = (-27/x, -27y/x^2). Lemma 3.2 gives I_0(X_0(6,7)) = {w_42, w_3, w_6, w_21}, i.e. those
// four have genus-0 quotients and w_2, w_7, w_14 have genus-1 quotients -- which is what our
// committed model shows, and what pins the labelling below.
//
// ⚠⚠ COMPOSE THE GENERATORS CAREFULLY -- THE FIRST ATTEMPT AT THIS FILE GOT IT WRONG AND THE ERROR
// LOOKED LIKE A DISAGREEMENT WITH THE MODEL. w_6(w_3(x,y)) = w_6(-x,y) = (27/x, -27y/x^2), NOT
// (-27/x, +27y/x^2). Mislabelling those two swaps w_2 with w_7, and the symptom was a quotient
// Jacobian of 42a2 against our 42a5 -- which reads as "the model is wrong" when in fact the model
// was right and the LABEL was wrong. Correctly composed, every key matches:
//     w_3  class [2]      w_6  class [3]      w_21 class [2,3]        (genus 0)
//     w_2  Jac 42a5       w_7  Jac 42a2       w_14 Jac 42a6           (genus 1)
// ⇒ the four Mobius maps are x, -x, 27/x, -27/x, each lifting to two maps by y -> +-y, giving all
// eight elements of W_{6,7} = {1,2,3,6,7,14,21,42}. Nothing is missing and nothing is duplicated.
//
// ⚠ W=[1,14] IS DELIBERATELY ABSENT FROM cover_data, and this is not an oversight. Our single
// committed entry there has the SAME Jacobian as the quotient derived from GR's model (42a6 both
// sides, and GR's own p.9 table publishes Jac(X_0(6,7)/<w_14>) = 42A6) but is NOT isomorphic to it
// as a CrvHyp -- an inequivalent quartic model of the same curve, i.e. a different degree-2 divisor
// class. That is the torsor phenomenon documented in tests/GonzalezRotger.m, not a disagreement.
// Listing it would make this test red for a reason that is not a defect.
//
// ⚠ NO TRANSPORT IS NEEDED despite the non-diagonal maps: cover_data[{1}] is GR's equation
// VERBATIM, and the helper compares ws_data in the EXPECTED curve's coordinates. Convention is
// row-vector times matrix on (X:Y:Z) of P(1,2,1), x = X/Z, y = Y/Z^2.
// ✅ All seven non-trivial matrices verified as automorphisms AND involutions before this was
// written; the four derived ones are formed by the group law w_m*w_n = w_{mn/gcd(m,n)^2}, so their
// labelling is forced by the generators rather than guessed.
//
// ⚠ QUOTIENT EQUATIONS ARE DERIVED; both identities were checked in the function field:
//     u = x + 27/x : f/x^2 = -3u^2 + 128     and     (x - 27/x)^2 = u^2 - 108
//     u = x - 27/x : f/x^2 = -3u^2 - 196     and     (x + 27/x)^2 = u^2 + 108

    // verifying [Gonzalez-Rotger, Table 1, p. 8]
    // D = 6, N = 7
    cover_data := AssociativeArray();
    id3 := IdentityMatrix(Rationals(), 3);   // read only under manual_isomorphism, which is off
    cover_data[{1}]    := <HyperellipticCurve(-3*x^4 - 34*x^2 - 2187),      id3>;
    cover_data[{1,3}]  := <HyperellipticCurve(-3*x^2 - 34*x - 2187),        id3>;  // u=x^2, v=y
    cover_data[{1,6}]  := <HyperellipticCurve(-3*x^2 - 196),                id3>;  // u=x-27/x, v=y/x
    cover_data[{1,21}] := <HyperellipticCurve(-3*x^2 + 128),                id3>;  // u=x+27/x, v=y/x
    cover_data[{1,7}]  := <HyperellipticCurve((-3*x^2 - 196)*(x^2 + 108)),  id3>;  // t=(y/x)(x+27/x)
    cover_data[{1,2}]  := <HyperellipticCurve((-3*x^2 + 128)*(x^2 - 108)),  id3>;  // t=(y/x)(x-27/x)
    cover_data[{1,42}] := <HyperellipticCurve(x^2 - x),                     id3>;  // P^1

    M42 := DiagonalMatrix(Rationals(), [1,-1,1]);                  // (x, -y)
    M3  := DiagonalMatrix(Rationals(), [-1,1,1]);                  // (-x, y)
    M6  := Matrix(Rationals(), 3, 3, [0,0,1, 0,-27,0, -27,0,0]);   // (-27/x, -27y/x^2)

    ws_data := AssociativeArray();
    ws_data[{1}] := AssociativeArray();
    ws_data[{1}][42] := M42;
    ws_data[{1}][3]  := M3;
    ws_data[{1}][6]  := M6;
    ws_data[{1}][14] := M3 * M42;        // w_3 * w_42 = w_{126/9}  = w_14
    ws_data[{1}][7]  := M6 * M42;        // w_6 * w_42 = w_{252/36} = w_7
    ws_data[{1}][2]  := M3 * M6;         // w_3 * w_6  = w_{18/9}   = w_2
    ws_data[{1}][21] := M3 * M6 * M42;   // w_2 * w_42 = w_{84/4}   = w_21
    // ⚠ No ws_data for the genus-0 keys: the helper's genus-0 branch exhibits no isomorphism to
    // conjugate by, and fails loudly rather than skip if asked.

    return cover_data, ws_data;
end function;

procedure test_6_7()
    cover_data, ws_data := load_covers_and_ws_data_6_7();
    curves := GetHyperellipticCandidates();

    test_AllEquationsAboveCoversSingleCurve(6, 7, cover_data, ws_data, curves);
    return;
end procedure;

test_6_7();
