import "tests/BorcherdsProducts.m" : test_AllEquationsAboveCoversSingleCurve;

// tests/X0_95_1.m -- RE-DERIVATION test for X_0^95(1).
//
// models_95_1.m is the first model of this base: its Borcherds search needed PR #63 (the second
// rung of the ladder had never finished before it).  What is checked:
//   [1] EXTERNAL: the stored W={1} curve against Guo-Yang's published equation (degree 16,
//       genus 7) -- in tests/GuoYangEquations.m, with an exact IsIsomorphic.
//   [2] RE-DERIVATION: AllEquationsAboveCovers is re-run and every cover below is compared.
//   [3] INVOLUTION: Guo-Yang's w_5 on the top curve, transported into our coordinates.
//
// Guo-Yang publish, in their coordinates, w_5(x,y) = (-1/x, y/x^8), a Mobius involution, so its
// matrix does not carry over as written (the argument of X0_69_1.m needs a diagonal change of
// coordinates, which this is not).  As in X0_55_1.m it was TRANSPORTED: psi := IsIsomorphic(our
// stored curve, Guo-Yang's) from the two equations alone, and the map recorded below is
// psi^-1 . w_5 . psi, which came out linear in the weighted coordinates (x : y : z):
//     (x, y, z) -> (-x - z, y, 2x + z)
// Its square is diag(-1, 1, -1) = the identity on the weighted ambient (weights 1, 8, 1), and it
// maps our curve to itself.  The other two involutions of the base are not transcribed in this
// repository yet and are not checked here.
//
// SOURCE for the equation and w_5: Compositio Math. 153 (2017) 1-40, Table A.1 "Equations of
// level one (continued)", printed page 35 (see tests/GuoYangEquations.m for the transcription).
//
// COST: the Borcherds forms take about 23 minutes on a Mac; the whole re-derivation about twice
// that on a loaded machine.

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
