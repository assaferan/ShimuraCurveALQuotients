import "tests/BorcherdsProducts.m" : test_AllEquationsAboveCoversSingleCurve;

function load_covers_and_ws_data_15_1()
    _<s> := PolynomialRing(Rationals());

// ✅ INVOLUTIONS CHECKED (2026-09-08), DERIVED from Guo-Yang's Example 32 rather than transcribed:
// that example gives no explicit involution list, but its construction forces one. They put s and
// y on X/w_3 (equation 3y^2 + (s+243)(s+3) = 0) and x on X/w_15 (via s = x^2, so x is the
// coordinate downstairs there). Hence w_3 must be the deck transformation of X -> X/w_3, which
// fixes s = x^2 and y and so is x -> -x; w_15 must fix x and negate y; and w_5 = w_3*w_15.
//     w_3(x,y) = (-x, y)    w_5(x,y) = (-x,-y)    w_15(x,y) = (x,-y)
// ⚠ NO TRANSPORT IS NEEDED HERE, unusually: cover_data below is already Guo-Yang's equation
// VERBATIM (3y^2 + (x^2+3)(x^2+243) = 0, their Example 32, agreeing with Jordan Prop. 3.2.1), so
// their coordinates ARE the expected curve's coordinates.
// ⚠ WHAT PINS THE LABELLING, since w_3 and w_15 both have genus-0 quotients and genus alone cannot
// separate them: it is WHICH quotient carries WHICH coordinate. X/w_3 is the one whose equation is
// in (s,y) -- and cover_data[{1,3}] is exactly -1/3*(s+243)*(s+3) -- while X/w_15 is the one with
// x^2 = s, and cover_data[{1,15}] is exactly s. Both already matched Guo-Yang before this change.
// Guo-Yang's own genus column (X/w_3 = 0, X/w_5 = 1, X/w_15 = 0) agrees with our model's keys.

    // verifying [Guo-Yang, Example 32, p. 22-24]
    // D = 15
    cover_data := AssociativeArray();
    cover_data[{1}] := <HyperellipticCurve(-1/3*(s^2+243)*(s^2+3)), DiagonalMatrix([9, 4*27, 1])>;
    cover_data[{1,3}] := <HyperellipticCurve(-1/3*(s+243)*(s+3)), DiagonalMatrix([-81/2, 4*27, 1])>;
    cover_data[{1,15}] := <HyperellipticCurve(s), DiagonalMatrix([-1/2, 1, 1]) >;

    ws_data := AssociativeArray();
    ws_data[{1}] := AssociativeArray();
    ws_data[{1}][3]  := DiagonalMatrix([-1, 1, 1]);
    ws_data[{1}][5]  := DiagonalMatrix([-1,-1, 1]);
    ws_data[{1}][15] := DiagonalMatrix([ 1,-1, 1]);

    return cover_data, ws_data;
end function;

procedure test_15_1()
    cover_data, ws_data := load_covers_and_ws_data_15_1();
    curves := GetHyperellipticCandidates();
    
    test_AllEquationsAboveCoversSingleCurve(15, 1, cover_data, ws_data, curves);
    return;
end procedure;

test_15_1();