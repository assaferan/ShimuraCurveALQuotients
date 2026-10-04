import "tests/BorcherdsProducts.m" : test_AllEquationsAboveCoversSingleCurve;

// tests/_offline/X0_95_1.m -- re-derivation test for X_0^95(1).
//
// Checks, against models_95_1.m:
//   [1] the stored W={1} curve is isomorphic to Guo-Yang's published equation (degree 16,
//       genus 7), in tests/GuoYangEquations.m;
//   [2] AllEquationsAboveCovers re-derives every cover and each agrees with the stored one;
//   [3] Guo-Yang's w_5 is an involution of the top curve, in our coordinates.
//
// Guo-Yang give w_5(x,y) = (-1/x, y/x^8) in their coordinates, a Mobius involution, so its matrix
// does not carry over as written.  As in X0_55_1.m it is transported: with psi the isomorphism
// from our stored curve to Guo-Yang's, obtained from the two equations alone, the map below is
// psi^-1 . w_5 . psi, which is linear in the weighted coordinates (x : y : z):
//     (x, y, z) -> (-x - z, y, 2x + z)
// Its square is diag(-1, 1, -1), the identity on the weighted ambient (weights 1, 8, 1), and it
// maps our curve to itself.  The other two involutions of the base are not transcribed in this
// repository and are not checked here.
//
// SOURCE for the equation and w_5: Compositio Math. 153 (2017) 1-40, Table A.1 "Equations of
// level one (continued)", printed page 35 (see tests/GuoYangEquations.m for the transcription).
//
// COST: 28 min on a Mac (the Borcherds forms are 23 min of it), so it is run by hand:
//     NORMALIZ_BIN=... magma -b filename:=tests/_offline/X0_95_1.m run_tests.m < /dev/null
// The stored models are checked in CI by tests/ModelChecks.m, including the L-polynomial of
// every quotient against the trace formula.

function load_covers_and_ws_data_95_1()
    _<s> := PolynomialRing(Rationals());

    cover_data := AssociativeArray();
    cover_data[{1,5}]  := <HyperellipticCurve(Polynomial(Rationals(),
        [ -7/3125, 132/15625, -1144/78125, 5766/390625, -6/625, 1628/390625, -467/390625,
          82/390625, -7/390625 ])), DiagonalMatrix([1,1,1])>;                               // genus 3
    cover_data[{1,95}] := <HyperellipticCurve(Polynomial(Rationals(),
        [ 1/5, -2/25, 1/25 ])), DiagonalMatrix([1,1,1])>;                                   // genus 0
    cover_data[{1}]    := <HyperellipticCurve(Polynomial(Rationals(),
        [ -7/390625, -138/390625, -1263/390625, -7086/390625, -27132/390625, -74764/390625,
          -152507/390625, -233628/390625, -270311/390625, -236122/390625, -154973/390625,
          -15192/78125, -27897/390625, -7954/390625, -1898/390625, -372/390625,
          -43/390625 ])), DiagonalMatrix([1,1,1])>;                                         // genus 7
    cover_data[{1,19}] := <HyperellipticCurve(Polynomial(Rationals(),
        [ -7/15625, 146/78125, -1443/390625, 8714/1953125, -36002/9765625, 21406/9765625,
          -9341/9765625, 2972/9765625, -666/9765625, 96/9765625, -7/9765625 ])),
        DiagonalMatrix([1,1,1])>;                                                           // genus 4

    // Guo-Yang's w_5, transported (see the header), as the harness applies it: (x, y, z) * M.
    ws_data := AssociativeArray();
    ws_data[{1}] := AssociativeArray();
    ws_data[{1}][5] := Matrix(3, 3, [ -1, 0, 2,   0, 1, 0,   -1, 0, 1 ]);
    return cover_data, ws_data;
end function;

procedure test_95_1()
    cover_data, ws_data := load_covers_and_ws_data_95_1();
    curves := GetHyperellipticCandidates();
    test_AllEquationsAboveCoversSingleCurve(95, 1, cover_data, ws_data, curves);
    return;
end procedure;

printf "testing equations of covers of X0*(95;1)...";
test_95_1();
printf " ok\n";
