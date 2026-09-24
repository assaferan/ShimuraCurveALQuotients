import "tests/BorcherdsProducts.m" : test_AllEquationsAboveCoversSingleCurve;

// tests/X0_22_5.m -- RE-DERIVATION test for X_0^22(5).
//
// Re-runs AllEquationsAboveCovers and compares what the pipeline actually produces, so passing IS
// reproduction -- the stronger claim than GuoYangEquations.m's stored-model comparison.
//
// ⚠ THIS BASE INCLUDES ITS `CRV` (paired) ENTRY, which was impossible before 2026-09-07. The
// helper now CONSTRUCTS the isomorphism for CRV pairs (tests/_crviso.m) instead of calling
// IsIsomorphic, which HANGS on them -- the genus-3 pair at 14_3 ran >50 min, the genus-5 one at
// 26_3 >1 h. The ambient weights are DERIVED (y's weight is half its own equation's degree), which
// tests/CRVStructure.m verifies against every stored entry.
//
// Keys with several stored entries are omitted: they cannot be matched unambiguously.
// The matrix in each cover_data value is a placeholder -- with manual_isomorphism false (the
// default) the helper never reads it.
//
// ✅ INVOLUTIONS CHECKED (2026-09-08), and the W={1} FULL CURVE is now compared -- both became
// possible only when EquationsByRebase landed. Before that the pipeline produced NOTHING for
// W={1}, {1,2}, {1,5}, {1,11}: it assembles a genus-g curve as a fibre product needing degree
// exactly g+1 over a shared base, and at this base the degrees produced were 1,2,4,6,7,8 with the
// genus-2 quotients needing 3. Changing the Hauptmodul on the star base fixes the degree profile.
//
// The expected curve below is GUO-YANG'S OWN published equation, so their involutions apply
// VERBATIM with no transport:
//     w_2(x,y) = (1/x, y/x^6)   w_5(x,y) = (-1/x, -y/x^6)   w_110(x,y) = (x, -y)
// On P(1,6,1) those are (x:y:z) -> (z:y:x), (-z:-y:x) and (x:-y:z) respectively.
// ⚠ NOT circular: the involutions and the curve are Guo-Yang's, and what the harness checks is
// that the PIPELINE's involution labelled w_m matches theirs, so a labelling error is detectable.

function load_covers_and_ws_data_22_5()
    _<s> := PolynomialRing(Rationals());


    cover_data := AssociativeArray();
    cover_data[{1}] := <HyperellipticCurve(Polynomial(Rationals(), [ -11, 0, -80, 0, -240, 0, -362, 0, -240, 0, -80, 0, -11 ])), DiagonalMatrix([1,1,1])>;   // genus 5 -- Guo-Yang's published equation
    cover_data[{1,10}] := <HyperellipticCurve(Polynomial(Rationals(), [ -1024/625, -4096/625, -6803/625, -6073/625, -3147/625, -951/625, -157/625, -11/625 ])), DiagonalMatrix([1,1,1])>;   // genus 3
    cover_data[{1,55}] := <HyperellipticCurve(Polynomial(Rationals(), [ -11/390625, 0, 6/78125, 0, 37/390625, 0, 56/390625, 0, 16/78125 ])), DiagonalMatrix([1,1,1])>;   // genus 3
    cover_data[{1,2,11,22}] := <HyperellipticCurve(Polynomial(Rationals(), [ -4096/625, 20044/625, -36799/625, 6008/125, -368/25 ])), DiagonalMatrix([1,1,1])>;   // genus 1
    cover_data[{1,2,5,10}] := <HyperellipticCurve(Polynomial(Rationals(), [ -1024/625, 951/125, -8267/625, 6376/625, -368/125 ])), DiagonalMatrix([1,1,1])>;   // genus 1
    cover_data[{1,5,22,110}] := <HyperellipticCurve(Polynomial(Rationals(), [ 0, -16, 16 ])), DiagonalMatrix([1,1,1])>;   // genus 0
    cover_data[{1,5,11,55}] := <HyperellipticCurve(Polynomial(Rationals(), [ 0, 1024/625, -3731/625, 4536/625, -368/125 ])), DiagonalMatrix([1,1,1])>;   // genus 1
    cover_data[{1,10,22,55}] := <HyperellipticCurve(Polynomial(Rationals(), [ 0, 65536/50625, -77248/10125, 909488/50625, -1069424/50625, 13952/1125, -5888/2025 ])), DiagonalMatrix([1,1,1])>;   // genus 2
    cover_data[{1,2,55,110}] := <HyperellipticCurve(Polynomial(Rationals(), [ 1, -9/4, 5/4 ])), DiagonalMatrix([1,1,1])>;   // genus 0
    cover_data[{1,10,11,110}] := <HyperellipticCurve(Polynomial(Rationals(), [ 0, -4, 5 ])), DiagonalMatrix([1,1,1])>;   // genus 0
    cover_data[{1,22}] := <HyperellipticCurve(Polynomial(Rationals(), [ -4096/625, 0, -732/125, 0, -1243/625, 0, -38/125, 0, -11/625 ])), DiagonalMatrix([1,1,1])>;   // genus 3

    ws_data := AssociativeArray();
    ws_data[{1}] := AssociativeArray();
    ws_data[{1}][2]   := Matrix(3,3,[ 0,0,1,  0, 1,0,  1,0,0 ]);
    ws_data[{1}][5]   := Matrix(3,3,[ 0,0,1,  0,-1,0, -1,0,0 ]);
    ws_data[{1}][110] := Matrix(3,3,[ 1,0,0,  0,-1,0,  0,0,1 ]);
    return cover_data, ws_data;
end function;

procedure test_22_5()
    cover_data, ws_data := load_covers_and_ws_data_22_5();
    curves := GetHyperellipticCandidates();
    test_AllEquationsAboveCoversSingleCurve(22, 5, cover_data, ws_data, curves);
    return;
end procedure;

test_22_5();
