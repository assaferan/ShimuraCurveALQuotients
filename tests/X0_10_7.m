import "tests/BorcherdsProducts.m" : test_AllEquationsAboveCoversSingleCurve;

function load_covers_and_ws_data_10_7()
    _<x> := PolynomialRing(Rationals());

// EXTERNAL SOURCE: Gonzalez-Rotger, "Non-elliptic Shimura curves of genus one",
// J. Math. Soc. Japan 58 (2006) 927-948; Table 1 p.8 and the involution table directly below it.
// Built on GONZALEZ-ROTGER, not Guo-Yang -- see tests/X0_14_1.m for the coverage gap this closes.
//
// ⚠⚠ THE (10,7) COLUMN OF THE p.8 INVOLUTION TABLE CONTAINS TWO ERRORS, BOTH CORRECTED HERE, BOTH
// CONFIRMED BY COMPUTATION. Read this before "restoring" either to what the paper prints.
//
//   (1) THE LABEL. The table prints `w_15`. But 15 does not divide D*N = 70, so w_15 is not an
//       Atkin-Lehner involution of X_0(10,7) at all -- W_{10,7} = {w_m : m | 70} =
//       {1,2,5,7,10,14,35,70}. Lemma 3.2 (p.7) gives I_0(X_0(10,7)) = {w_70, w_5, w_10, w_35},
//       so the intended label is w_5 or w_35. It is **w_5**: the quotient of GR's curve by
//       (x,y) -> (-1/x, -y/x^2) is the conic -27u^2 - 40u - 48, whose class is [2], and our
//       committed W=[1,5] has class [2] in all three of its entries, while W=[1,35] is [2,5] and
//       W=[1,10] is [5]. The class discriminates, and it picks w_5.
//
//   (2) THE MOBIUS MAP. The table prints w_10(x,y) = ((2x-1)/(x-2), 5y/(x-2)^2). With `2x-1` this
//       is not even a self-map of the curve: (x-2)^4 * f((2x-1)/(x-2)) / f(x) is not constant, it
//       is a ratio of quartics. It is also not an involution -- x' - 2 = 3/(x-2) there, so
//       y'' = 25y/9 rather than y. The sign is wrong: with **(2x+1)/(x-2)** we get x' - 2 =
//       5/(x-2), hence y'' = y, and the published y-part `5y/(x-2)^2` is then exactly right.
//       ⇒ CONFIRMED INDEPENDENTLY: IsGL2Equivalent(f, f, 4) returns exactly four self-equivalences
//       of GR's quartic -- x, -1/x, (2x+1)/(x-2) and their composite -- and (2x-1)/(x-2) is not
//       among them. Four Mobius maps, each lifting to two maps on the curve, is the full group of
//       order 8 = #W_{10,7}, so nothing is missing.
//
// ✅ ALL SEVEN NON-TRIVIAL INVOLUTIONS VERIFIED as automorphisms of GR's curve AND as involutions,
// by IsIsomorphism and map equality, before this file was written. The three generators are written
// out; the rest are formed by the Atkin-Lehner group law w_m * w_n = w_{mn/gcd(m,n)^2}, so the
// labelling of the derived ones is forced by the labelling of the generators rather than guessed.
//
// ⚠ NO TRANSPORT IS NEEDED, despite the non-diagonal maps: cover_data[{1}] is GR's equation
// VERBATIM, and the helper compares ws_data in the EXPECTED curve's coordinates. Convention is
// row-vector times matrix on (X:Y:Z) of P(1,2,1), with x = X/Z and y = Y/Z^2.
//
// ⚠ THE QUOTIENT EQUATIONS ARE DERIVED, and were checked in the function field and then against the
// committed model:
//     w_5 : u = x - 1/x, v = y/x invariant, v^2 = f/x^2 = -27u^2 - 40u - 48   class [2]  = ours
//     w_10: u = (x^2+1)/(x-2), t = y/(x-2) invariant, t^2 = f/(x-2)^2
//                                              = -27u^2 - 40u - 20           class [5]  = ours
//     w_70: the x-line, P^1 over Q by their Lemma 2.1                        split      = ours

    // verifying [Gonzalez-Rotger, Table 1, p. 8]
    // D = 10, N = 7
    cover_data := AssociativeArray();
    // ⚠ The second component is read only when manual_isomorphism is set, which it is not here.
    id3 := IdentityMatrix(Rationals(), 3);
    cover_data[{1}]    := <HyperellipticCurve(-27*x^4 - 40*x^3 + 6*x^2 + 40*x - 27), id3>;
    cover_data[{1,5}]  := <HyperellipticCurve(-27*x^2 - 40*x - 48),                  id3>;
    cover_data[{1,10}] := <HyperellipticCurve(-27*x^2 - 40*x - 20),                  id3>;
    cover_data[{1,70}] := <HyperellipticCurve(x^2 - x),                              id3>;  // P^1

    // the three generators, in GR's own coordinates
    M70 := DiagonalMatrix(Rationals(), [1,-1,1]);                // (x,y) -> (x, -y)
    M5  := Matrix(Rationals(), 3, 3, [0,0,1, 0,-1,0, -1,0,0]);   // (x,y) -> (-1/x, -y/x^2)
    M10 := Matrix(Rationals(), 3, 3, [2,0,1, 0,5,0, 1,0,-2]);    // (x,y) -> ((2x+1)/(x-2), 5y/(x-2)^2)

    ws_data := AssociativeArray();
    ws_data[{1}] := AssociativeArray();
    ws_data[{1}][70] := M70;
    ws_data[{1}][5]  := M5;
    ws_data[{1}][10] := M10;
    ws_data[{1}][14] := M5 * M70;          // w_5  * w_70 = w_{350/25}  = w_14
    ws_data[{1}][7]  := M10 * M70;         // w_10 * w_70 = w_{700/100} = w_7
    ws_data[{1}][2]  := M5 * M10;          // w_5  * w_10 = w_{50/25}   = w_2
    ws_data[{1}][35] := M5 * M10 * M70;    // w_2  * w_70 = w_{140/4}   = w_35
    // ⚠ No ws_data for the genus-0 keys: the helper's genus-0 branch decides by conic class and
    // exhibits no isomorphism to conjugate by, and fails loudly rather than skip if asked.

    return cover_data, ws_data;
end function;

procedure test_10_7()
    cover_data, ws_data := load_covers_and_ws_data_10_7();
    curves := GetHyperellipticCandidates();

    test_AllEquationsAboveCoversSingleCurve(10, 7, cover_data, ws_data, curves);
    return;
end procedure;

test_10_7();
