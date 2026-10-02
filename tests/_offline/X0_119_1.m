import "tests/BorcherdsProducts.m" : test_AllEquationsAboveCoversSingleCurve;

// tests/_offline/X0_119_1.m -- re-derivation test for X_0^119(1).
//
// Checks, against models_119_1.m:
//   [1] the stored W={1} curve is isomorphic to Guo-Yang's published equation (degree 20,
//       genus 9), in tests/GuoYangEquations.m;
//   [2] AllEquationsAboveCovers re-derives every cover and each agrees with the stored one;
//   [3] Guo-Yang's w_7 is an involution of the top curve.
//
// Guo-Yang give w_7(x,y) = (-x, y) in their coordinates.  Our stored W={1} polynomial is even
// in x, as is Guo-Yang's, so the two models differ by a diagonal change of coordinates and a
// sign-only involution has the same matrix in both (the argument of X0_69_1.m).  The other two
// involutions of the base are not transcribed in this repository and are not checked here.
//
// SOURCE for the equation and w_7: Compositio Math. 153 (2017) 1-40, Table A.1 "Equations of
// level one (continued)", printed page 35 (see tests/GuoYangEquations.m for the transcription).
//
// COST: 58 min on a Mac, so it is run by hand:
//     NORMALIZ_BIN=... magma -b filename:=tests/_offline/X0_119_1.m run_tests.m < /dev/null
// The stored models are checked in CI by tests/ModelChecks.m, including the L-polynomial of
// every quotient against the trace formula.

function load_covers_and_ws_data_119_1()
    _<s> := PolynomialRing(Rationals());

    cover_data := AssociativeArray();
    cover_data[{1,7}] := <HyperellipticCurve(Polynomial(Rationals(),
        [ -1/2121269248, -533/363797676032, -74475/35652172251136, -368239/218369555038208,
          -13343/17826086125568, -9113/62391301439488, 5057/873478220152832,
          -57/31195650719744, -61/249565205757952, 75/873478220152832, -1/249565205757952 ])), DiagonalMatrix([1,1,1])>;   // genus 4
    cover_data[{1,17}] := <HyperellipticCurve(Polynomial(Rationals(),
        [ 0, -343/303038464, -533/151519232, -74475/14848884736, -368239/90949419008,
          -13343/7424442368, -9113/25985548288, 5057/363797676032, -57/12992774144,
          -61/103942193152, 75/363797676032, -1/103942193152 ])), DiagonalMatrix([1,1,1])>;   // genus 5
    cover_data[{1,119}] := <HyperellipticCurve(Polynomial(Rationals(),
        [ 0, 1/49 ])), DiagonalMatrix([1,1,1])>;   // genus 0
    cover_data[{1}] := <HyperellipticCurve(Polynomial(Rationals(),
        [ -1/2121269248, 0, -533/7424442368, 0, -74475/14848884736, 0, -368239/1856110592, 0,
          -653807/151519232, 0, -3125759/75759616, 0, 12141857/151519232, 0,
          -46941951/37879808, 0, -2461570027/303038464, 0, 21185643675/151519232, 0,
          -96889010407/303038464 ])), DiagonalMatrix([1,1,1])>;   // genus 9

    ws_data := AssociativeArray();
    ws_data[{1}] := AssociativeArray();
    ws_data[{1}][7] := DiagonalMatrix([-1, 1, 1]);
    return cover_data, ws_data;
end function;

procedure test_119_1()
    cover_data, ws_data := load_covers_and_ws_data_119_1();
    curves := GetHyperellipticCandidates();
    test_AllEquationsAboveCoversSingleCurve(119, 1, cover_data, ws_data, curves);
    return;
end procedure;

printf "testing equations of covers of X0*(119;1)...";
test_119_1();
printf " ok\n";
