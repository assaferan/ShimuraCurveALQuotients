import "tests/BorcherdsProducts.m" : test_AllEquationsAboveCoversSingleCurve;

// tests/_offline/X0_159_1.m -- re-derivation test for X_0^159(1).
//
// Checks, against models_159_1.m:
//   [1] the stored W={1} curve is isomorphic to Guo-Yang's published equation (degree 20,
//       genus 9), in tests/GuoYangEquations.m;
//   [2] AllEquationsAboveCovers re-derives every cover and each agrees with the stored one;
//   [3] Guo-Yang's w_3 is an involution of the top curve.
//
// Guo-Yang give w_3(x,y) = (-x, y) in their coordinates.  Our stored W={1} polynomial is even
// in x, as is Guo-Yang's, so the two models differ by a diagonal change of coordinates and a
// sign-only involution has the same matrix in both (the argument of X0_69_1.m).  The other two
// involutions of the base are not transcribed in this repository and are not checked.  The
// [1,53] quotient is stored as y^2 + h y = f.
//
// SOURCE for the equation and w_3: Compositio Math. 153 (2017) 1-40, Table A.1 "Equations of
// level one (continued)", printed page 35 (see tests/GuoYangEquations.m for the transcription).
//
// COST: 10395 s (2.9 h) on a Mac, the third rung's kernel dominating, so it is run by hand:
//     NORMALIZ_BIN=... magma -b filename:=tests/_offline/X0_159_1.m run_tests.m < /dev/null
// The stored models are checked in CI by tests/ModelChecks.m, including the L-polynomial of
// every quotient against the trace formula.

function load_covers_and_ws_data_159_1()
    _<s> := PolynomialRing(Rationals());

    cover_data := AssociativeArray();
    cover_data[{1,3}] := <HyperellipticCurve(Polynomial(Rationals(),
        [ -14348907/1048576, 320458923/524288, -2137987143/1048576, -19356675543/131072,
          271620026541/524288, 1448708697441/262144, -15903654120171/524288,
          -5314371204807/131072, -20887603477551/1048576, -2707681797621/524288,
          -847288609443/1048576 ])), DiagonalMatrix([1,1,1])>;   // genus 4
    cover_data[{1,53}] := <HyperellipticCurve(Polynomial(Rationals(),
        [ -5, 73, -437, 1221, -1361, -286, 1063, 571, -483, -46, 25, 3 ]), Polynomial(Rationals(),
        [ 1, 0, 1, 1, 1, 1 ])), DiagonalMatrix([1,1,1])>;   // genus 5, y^2 + h*y = f
    cover_data[{1,159}] := <HyperellipticCurve(Polynomial(Rationals(),
        [ 0, 9 ])), DiagonalMatrix([1,1,1])>;   // genus 0
    cover_data[{1}] := <HyperellipticCurve(Polynomial(Rationals(),
        [ -14348907/1048576, 0, 35606547/524288, 0, -26394903/1048576, 0, -26552367/131072, 0,
          41399181/524288, 0, 24534009/262144, 0, -29925531/524288, 0, -1111103/131072, 0,
          -485231/1048576, 0, -6989/524288, 0, -243/1048576 ])), DiagonalMatrix([1,1,1])>;   // genus 9

    ws_data := AssociativeArray();
    ws_data[{1}] := AssociativeArray();
    ws_data[{1}][3] := DiagonalMatrix([-1, 1, 1]);
    return cover_data, ws_data;
end function;

procedure test_159_1()
    cover_data, ws_data := load_covers_and_ws_data_159_1();
    curves := GetHyperellipticCandidates();
    test_AllEquationsAboveCoversSingleCurve(159, 1, cover_data, ws_data, curves);
    return;
end procedure;

printf "testing equations of covers of X0*(159;1)...";
test_159_1();
printf " ok\n";
