// ⚠ OFFLINE, ON COST GROUNDS ONLY (`run_tests.m` globs just `tests/*.m`). Re-deriving 39_2 takes
// over 17 minutes, which does not belong in a CI slot; the test itself is expected to pass.
// It became possible at all on 2026-09-07, when the coprime-to-level CM filter became OFF BY
// DEFAULT: before that 39_2 needed CMNONCOPRIME=1, and a plain X0_*.m test cannot set an env var.
//
import "tests/BorcherdsProducts.m" : test_AllEquationsAboveCoversSingleCurve;

// tests/X0_39_2.m -- RE-DERIVATION test for X_0^39(2).
//
// Unlike tests/GuoYangEquations.m, which reads the STORED model, this re-runs
// AllEquationsAboveCovers and compares what the pipeline actually produces. Passing therefore
// IS reproduction, which is the stronger claim (see data/models/PROVENANCE.md).
//
// The expected curves below are the committed model's hyperelliptic cover entries, which are
// themselves validated against Guo-Yang's published equations. CRV (paired) entries and empty
// keys are omitted: the helper compares hyperelliptic covers, and a key with several stored
// entries cannot be matched unambiguously.
//
// The second component of each cover_data value is unused here -- with manual_isomorphism false
// (the default) the helper calls IsIsomorphic, so the matrix is a placeholder. ws_data is left
// empty for the same reason: the helper skips involution checks for keys it does not find.

function load_covers_and_ws_data_39_2()
    _<s> := PolynomialRing(Rationals());

    cover_data := AssociativeArray();
    cover_data[{1,13}] := <HyperellipticCurve(Polynomial(Rationals(), [ -27/4096, 0, -1585/173056, 0, -641/346112, 0, -25/173056, 0, -3/692224 ])), DiagonalMatrix([1,1,1])>;   // genus 3
    cover_data[{1}] := <HyperellipticCurve(Polynomial(Rationals(), [ -63/169, 90/169, -261/169, 900/169, -1161/169, 360/169, -4239/169, -450/169, -360/13, 450/169, -4239/169, -360/169, -1161/169, -900/169, -261/169, -90/169, -63/169 ])), DiagonalMatrix([1,1,1])>;   // genus 7
    cover_data[{1,6,26,39}] := <HyperellipticCurve(Polynomial(Rationals(), [ 0, -81/2, 81/2 ])), DiagonalMatrix([1,1,1])>;   // genus 0
    cover_data[{1,3}] := <HyperellipticCurve(Polynomial(Rationals(), [ -243/8192, 0, -61623/1384448, 0, -8939/692224, 0, -1091/692224, 0, -127/1384448, 0, -3/1384448 ])), DiagonalMatrix([1,1,1])>;   // genus 4
    cover_data[{1,2}] := <HyperellipticCurve(Polynomial(Rationals(), [ -81/3328, 0, -549/21632, 0, 193/21632, 0, -7/5408, 0, 35/43264, 0, -3/21632 ])), DiagonalMatrix([1,1,1])>;   // genus 4
    cover_data[{1,2,13,26}] := <HyperellipticCurve(Polynomial(Rationals(), [ -81/3328, -4941/43264, 15633/86528, -5103/43264, 229635/692224, -177147/692224 ])), DiagonalMatrix([1,1,1])>;   // genus 2
    cover_data[{1,6}] := <HyperellipticCurve(Polynomial(Rationals(), [ 9/3328, 9/416, 63/2704, -441/2704, -495/1352, 9/338, 153/676, -45/169, -63/169 ])), DiagonalMatrix([1,1,1])>;   // genus 3
    cover_data[{1,3,13,39}] := <HyperellipticCurve(Polynomial(Rationals(), [ -9, 9 ])), DiagonalMatrix([1,1,1])>;   // genus 0
    cover_data[{1,2,3,6}] := <HyperellipticCurve(Polynomial(Rationals(), [ 0, -9477/512, -44469/512, 140697/1024, -45927/512, 2066715/8192, -1594323/8192 ])), DiagonalMatrix([1,1,1])>;   // genus 2
    cover_data[{1,26}] := <HyperellipticCurve(Polynomial(Rationals(), [ -81/3328, -405/1664, -24867/43264, 6723/5408, 17415/2704, 13203/2704, -7857/1352, -243/169, 6885/676, 729/169, -567/169 ])), DiagonalMatrix([1,1,1])>;   // genus 4
    cover_data[{1,6,13,78}] := <HyperellipticCurve(Polynomial(Rationals(), [ 9/3328, 333/21632, -405/86528, 729/86528, -19683/692224 ])), DiagonalMatrix([1,1,1])>;   // genus 1
    cover_data[{1,2,39,78}] := <HyperellipticCurve(Polynomial(Rationals(), [ 0, 9/2 ])), DiagonalMatrix([1,1,1])>;   // genus 0
    cover_data[{1,3,26,78}] := <HyperellipticCurve(Polynomial(Rationals(), [ 0, 81/6656, 2997/43264, -3645/173056, 6561/173056, -177147/1384448 ])), DiagonalMatrix([1,1,1])>;   // genus 2
    cover_data[{1,78}] := <HyperellipticCurve(Polynomial(Rationals(), [ 9/3328, 0, 37/10816, 0, -5/21632, 0, 1/10816, 0, -3/43264 ])), DiagonalMatrix([1,1,1])>;   // genus 3

    ws_data := AssociativeArray();
    return cover_data, ws_data;
end function;

procedure test_39_2()
    cover_data, ws_data := load_covers_and_ws_data_39_2();
    curves := GetHyperellipticCandidates();
    test_AllEquationsAboveCoversSingleCurve(39, 2, cover_data, ws_data, curves);
    return;
end procedure;

test_39_2();
