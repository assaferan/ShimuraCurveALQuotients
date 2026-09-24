import "tests/BorcherdsProducts.m" : test_AllEquationsAboveCoversSingleCurve;

// tests/X0_69_1.m -- RE-DERIVATION test for X_0^69(1).
//
// ⚠ WHY THIS EXISTS, and it is not "one more base". Until 2026-09-22 `69_1` was validated by
// NOTHING external. It is a Guo-Yang base with a PUBLISHED equation, yet it appeared neither in
// tests/GuoYangEquations.m nor in any X0_*.m test, so the only thing checking it was ModelChecks
// -- structural self-consistency, which cannot see a wrong curve.
//
// The cause was pure timing, and it is worth recording because it will recur. models_69_1.m landed
// 2026-09-14 in `4bfb859` (ScaleForSchofer: stop treating w_1 as an Atkin-Lehner involution, the
// fix that built 33_1 and 69_1). Both X0_* generation batches predate it -- 2025-11-19 and
// 2026-09-06/07 -- and the equation table was last extended before it too. A model that arrives
// after the sweep is never swept, and nothing in the repo asks "which committed models have no
// oracle?". ⇒ After any fix that BUILDS a base, check whether it is a published base and wire up
// both halves; a model with no oracle is the condition that let the 10_3 conic drift live a week.
//
// WHAT IS CHECKED:
//   [1] EXTERNAL: the stored W={1} curve against Guo-Yang's published
//           y^2 = -243x^8 + 1268x^6 - 666x^4 - 2268x^2 - 2187            [degree 8, genus 3]
//       That comparison lives in tests/GuoYangEquations.m (added the same day, exact
//       IsIsomorphic, two perturbation controls). This file adds the RE-DERIVATION on top.
//   [2] RE-DERIVATION: AllEquationsAboveCovers is re-run and every cover below is compared.
//   [3] INVOLUTIONS: Guo-Yang's own w_3 and w_69 on the top curve -- see the note on ws_data.
//
// ✅ NEGATIVE-CONTROLLED 2026-09-22, and the control is what makes [3] evidence. Swapping the two
// matrices (w_3 <-> w_69) makes this test FAIL:
//     "no isomorphism ... intertwines all 2 labelled Atkin-Lehner involution(s). 5 candidate
//      map(s) tried. The curves ARE isomorphic -- what fails is that no identification matches
//      the involution LABELLING."
// Note what that message establishes: 2 involutions were actually compared (not silently skipped),
// the helper searched FIVE maps in the Isom torsor before giving up rather than passing on a lucky
// phi, and the curves remained isomorphic throughout -- so the check discriminates the LABELLING
// specifically, which is the only thing that makes these quotient models rather than curves.
// Re-run the control by swapping the two DiagonalMatrix lines below; it must go red.
//
// COST 117.5 s, so this file belongs in tests/ and NOT tests/_offline/ -- it is checked on every
// push. ⚠ That makes 69_1 the only one of 69_1/87_1/39_2/111_1/93_1 that CI sees at all; the
// other four are offline and run only when someone remembers them.
//
// SOURCE for the equation and the involutions: Compositio Math. 153 (2017) 1-40, Table A.1
// "Equations of level one (continued)", printed page 34. Read three ways (journal page visually,
// journal PDF text layer, arXiv v1 TeX), all agreeing.

function load_covers_and_ws_data_69_1()
    _<s> := PolynomialRing(Rationals());

    // The committed model's hyperelliptic cover entries, which are themselves now validated
    // against Guo-Yang's published equation (see [1] above).
    cover_data := AssociativeArray();
    cover_data[{1,3}]  := <HyperellipticCurve(Polynomial(Rationals(),
        [ -4/6561, 1/8748, 35/139968, -37/839808, -1/27648 ])), DiagonalMatrix([1,1,1])>;        // genus 1
    cover_data[{1,23}] := <HyperellipticCurve(Polynomial(Rationals(),
        [ -4/81, 19/324, 19/1728, -247/10368, 53/82944, 3/1024 ])), DiagonalMatrix([1,1,1])>;    // genus 2
    cover_data[{1,69}] := <HyperellipticCurve(Polynomial(Rationals(),
        [ 9, -9 ])), DiagonalMatrix([1,1,1])>;                                                   // genus 0
    cover_data[{1}]    := <HyperellipticCurve(Polynomial(Rationals(),
        [ -1/3072, 0, -7/186624, 0, -37/30233088, 0, 317/1224440064, 0, -1/181398528 ])),
        DiagonalMatrix([1,1,1])>;                                                                // genus 3

    // ATKIN-LEHNER INVOLUTIONS, transcribed from Guo-Yang (journal, printed page 34):
    //     w_3 (x,y) = (-x,  y)        w_69 (x,y) = ( x, -y)
    // ⚠ WHY THESE CARRY OVER UNCHANGED into our coordinates -- the same argument X0_51_1.m makes,
    // and it has to be checked per base rather than assumed. BOTH our stored W={1} polynomial and
    // Guo-Yang's are EVEN in x (every odd coefficient of the stored entry above is literally 0),
    // so the two models differ by a DIAGONAL change x -> ax, y -> by. Diagonal scalings commute
    // with sign changes, so a sign-only involution has the same matrix in both coordinate systems.
    // This would NOT hold for a Mobius involution such as 111_1's w_37 = (-1/x, y/x^8); those are
    // still expressible as a matrix on the weighted ambient (X0_55_1.m does exactly that for
    // w_5 = (-1/x, y/x^4)), but the matrix must then be CONJUGATED through the isomorphism rather
    // than copied across.
    // The ambient of a HyperellipticCurve has coordinates (x, y, z); z is unaffected.
    ws_data := AssociativeArray();
    ws_data[{1}] := AssociativeArray();
    ws_data[{1}][3]  := DiagonalMatrix([-1, 1, 1]);
    ws_data[{1}][69] := DiagonalMatrix([ 1,-1, 1]);
    return cover_data, ws_data;
end function;

procedure test_69_1()
    cover_data, ws_data := load_covers_and_ws_data_69_1();
    curves := GetHyperellipticCandidates();
    test_AllEquationsAboveCoversSingleCurve(69, 1, cover_data, ws_data, curves);
    return;
end procedure;

printf "testing equations of covers of X0*(69;1)...";
test_69_1();
printf " ok\n";
