// ✅ DIAGNOSED AND FIXED 2026-09-09. This test used to FAIL with the cause unknown; it was a
// DROPPED h-TERM in the generated `cover_data`, not a defect in the pipeline or the model.
//
// models_87_1.m stores its [1,29] entry as `<3, f, h>` with `h = x^3 + x^2 + 1`, meaning the curve
// is `y^2 + h*y = f`. The generator (scratch, `genx0.py`) emitted only `f`, so the test compared
// `y^2 = f` against what the pipeline produces -- a DIFFERENT curve, of the same genus, so the
// genus check could not catch it. Nine entries across seven model files carry a nonzero `h`; this
// was the only generated test affected.
//
// ⚠ THE PREVIOUS HEADER'S DIAGNOSIS WAS SOUND AND POINTED THE RIGHT WAY: it established that the
// model regenerates exactly (`ModelRegen`: "87_1: OK") and that `GuoYangEquations.m` independently
// validates the `W={1}` curve, and concluded the defect was "in THIS generated test, not in the
// data it compares against". That was correct. What was missing was only the mechanism.
//
// ⚠ IT REMAINS OFFLINE because it is SLOW (3945 s = 66 min), not because it is broken.
// `run_tests.m` globs only `tests/*.m`, so `_offline` cannot reach CI.
//
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
    // ⚠ THE h-TERM. models_87_1.m stores this entry as <3, f, h> with h = x^3 + x^2 + 1, i.e. the
    // curve is y^2 + h*y = f, NOT y^2 = f. The generator dropped h, so this line used to compare a
    // DIFFERENT curve of the same genus -- which is why this test failed. See the header.
    cover_data[{1,29}] := <HyperellipticCurve(Polynomial(Rationals(), [ -5, 14, 23, -81, -36, 93, 70, 18, 3 ]), Polynomial(Rationals(), [ 1, 0, 1, 1 ])), DiagonalMatrix([1,1,1])>;   // genus 3
    cover_data[{1,3}] := <HyperellipticCurve(Polynomial(Rationals(), [ -129140163/3444736, 1190959281/1722368, -15635525661/3444736, 10581521751/861184, -34231709133/3444736, -8451506223/1722368, -10460353203/3444736 ])), DiagonalMatrix([1,1,1])>;   // genus 2
    cover_data[{1,87}] := <HyperellipticCurve(Polynomial(Rationals(), [ 0, -27 ])), DiagonalMatrix([1,1,1])>;   // genus 0
    cover_data[{1}] := <HyperellipticCurve(Polynomial(Rationals(), [ -129140163/3444736, 0, -44109603/1722368, 0, -21447909/3444736, 0, -537597/861184, 0, -64413/3444736, 0, 589/1722368, 0, -27/3444736 ])), DiagonalMatrix([1,1,1])>;   // genus 5

    // ATKIN-LEHNER INVOLUTIONS, transcribed from Guo-Yang (journal, Table A.1, printed page 34):
    //     w_3 (x,y) = (-x,  y)        w_87 (x,y) = ( x, -y)
    // ✅ ADDED 2026-09-23. Until now this file's ws_data was EMPTY -- boilerplate inherited from
    // the 2026-09-07 generated batch ("the helper skips involution checks for keys it does not
    // find"), i.e. undone work rather than a decision. A pass without it says the CURVES are
    // re-derived, not that their LABELLING is right, and the labelling is what makes these
    // quotient models rather than merely curves.
    //
    // ⚠ WHY THESE CARRY OVER UNCHANGED, checked and not assumed: the stored W={1} polynomial above
    // is EVEN in x (every odd coefficient is literally 0), as is Guo-Yang's published
    // y^2 = -(x^6-7x^4+43x^2+27)(243x^6+523x^4+369x^2+81). So the two differ by a DIAGONAL change
    // x -> ax, y -> by, and diagonal scalings commute with sign changes, leaving a sign-only
    // involution with the same matrix on both sides. Same argument as X0_51_1.m and X0_69_1.m.
    // ⚠ VERIFIED AS AUTOMORPHISMS of the stored curve before this file was run at all (the run is
    // hours): both matrices below are automorphisms AND involutions of cover_data[{1}], while the
    // x <-> z swap Matrix(3,3,[0,0,1, 0,1,0, -1,0,0]) is NOT -- which is the discriminating fact,
    // since at 111_1 the situation is exactly reversed.
    ws_data := AssociativeArray();
    ws_data[{1}] := AssociativeArray();
    ws_data[{1}][3]  := DiagonalMatrix([-1, 1, 1]);
    ws_data[{1}][87] := DiagonalMatrix([ 1,-1, 1]);
    return cover_data, ws_data;
end function;

procedure test_87_1()
    cover_data, ws_data := load_covers_and_ws_data_87_1();
    curves := GetHyperellipticCandidates();
    test_AllEquationsAboveCoversSingleCurve(87, 1, cover_data, ws_data, curves);
    return;
end procedure;

test_87_1();
