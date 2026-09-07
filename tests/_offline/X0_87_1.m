// ⚠⚠ THIS TEST CURRENTLY FAILS, AND THE CAUSE IS NOT YET KNOWN. It is kept OFFLINE
// (`run_tests.m` globs only `tests/*.m`) so it cannot reach CI, and committed so the open question
// is not lost.
//
// WHAT IS KNOWN, 2026-09-07:
//   * the MODEL is fine -- `tests/_offline/ModelRegen.m` reports `87_1: OK (4 covers compared,
//     0 CRV skipped)`, i.e. `models_87_1.m` regenerates from the pipeline exactly;
//   * `tests/GuoYangEquations.m` independently checks `87_1`'s `W={1}` genus-5 curve against
//     Guo-Yang's published equation and passes;
//   * so the defect is in THIS generated test, not in the data it compares against.
//   * it is also slow: 3945 s (66 min), which is why it would not belong in CI even once fixed.
//
// The generator that produced it (scratch, `genx0.py`) emits the committed model's hyperelliptic
// cover entries as `cover_data` and lets the helper call `IsIsomorphic`. For the other four bases
// generated the same way -- `51_1`, `55_1`, `57_1`, `14_5` -- that works and they are in CI.
// ⇒ NEXT STEP when picked up: run it with `verbose:=1` and read the coverage line the helper now
// prints ("N curve comparison(s) ... M/K expected covers matched"). That distinguishes "an
// expected cover was never produced" from "a produced cover failed its isomorphism", which are
// different problems.

import "tests/BorcherdsProducts.m" : test_AllEquationsAboveCoversSingleCurve;

// tests/X0_87_1.m -- RE-DERIVATION test for X_0^87(1).
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

function load_covers_and_ws_data_87_1()
    _<s> := PolynomialRing(Rationals());

    cover_data := AssociativeArray();
    cover_data[{1,29}] := <HyperellipticCurve(Polynomial(Rationals(), [ -5, 14, 23, -81, -36, 93, 70, 18, 3 ])), DiagonalMatrix([1,1,1])>;   // genus 3
    cover_data[{1,3}] := <HyperellipticCurve(Polynomial(Rationals(), [ -129140163/3444736, 1190959281/1722368, -15635525661/3444736, 10581521751/861184, -34231709133/3444736, -8451506223/1722368, -10460353203/3444736 ])), DiagonalMatrix([1,1,1])>;   // genus 2
    cover_data[{1,87}] := <HyperellipticCurve(Polynomial(Rationals(), [ 0, -27 ])), DiagonalMatrix([1,1,1])>;   // genus 0
    cover_data[{1}] := <HyperellipticCurve(Polynomial(Rationals(), [ -129140163/3444736, 0, -44109603/1722368, 0, -21447909/3444736, 0, -537597/861184, 0, -64413/3444736, 0, 589/1722368, 0, -27/3444736 ])), DiagonalMatrix([1,1,1])>;   // genus 5

    ws_data := AssociativeArray();
    return cover_data, ws_data;
end function;

procedure test_87_1()
    cover_data, ws_data := load_covers_and_ws_data_87_1();
    curves := GetHyperellipticCandidates();
    test_AllEquationsAboveCoversSingleCurve(87, 1, cover_data, ws_data, curves);
    return;
end procedure;

test_87_1();
