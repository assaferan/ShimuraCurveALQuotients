// ⚠ MOVED OUT OF CI 2026-09-06 -- IT WAS VERIFYING NOTHING THERE.
//
// This test passes LOCALLY with real comparisons, but in GitHub CI it made ZERO curve
// comparisons: `AllEquationsAboveCovers` produced the expected `W={1}` key with NO BASES, so the
// comparison loop never ran. It had been passing green on that basis, and cost 5016 s (84 min) per
// run to do it. The vacuity guard added to `test_AllEquationsAboveCoversSingleCurve` on
// 2026-09-06 is what exposed it.
//
// ROOT CAUSE: **CI never sets `NORMALIZ_BIN`.** `CLAUDE.md` is explicit that without it a fresh
// polytope solve fails SILENTLY -- "you get 'no solutions' rather than an error" -- so any cover
// needing a solve beyond the committed cache simply comes back empty. Locally the variable is set
// and the same cover is found.
//
// ⚠ MEASURED, so the scope is known: every OTHER X0_* job in CI reports full coverage
// (`X0_10_11` 1/1, `X0_26_1` 3/3, `X0_6_17` 1/1), so `10_19` is the only test affected. This is
// not a general CI collapse.
//
// ⇒ THE REAL FIX is to install Normaliz in CI and set `NORMALIZ_BIN`, after which this file should
// move back to `tests/`. Until then it lives here, where it runs meaningfully.
//
import "tests/BorcherdsProducts.m" : test_AllEquationsAboveCoversSingleCurve;

function load_covers_and_ws_data_10_19()
    P3<x,y,z,s> := WeightedProjectiveSpace(Rationals(), [1,3,1,1]);
     // D = 82
    cover_data := AssociativeArray();
    cover_data[{1}] := <Curve(P3, [y^2 + 8*x^6 - 57*x^4*s^2 + 40*x^2*s^4 - 16*s^6, z^2 - 5*x^2 + 32*s^2]), Matrix([[0,0,1,0],[0,1/8,0,0],[-1/8,0,0,-1/8],[-1/4,0,0,0]])>;

    ws_data := AssociativeArray();
    ws_data[{1}] := AssociativeArray();
    ws_data[{1}][2] := DiagonalMatrix([-1,1,1,1]);
    ws_data[{1}][5] := DiagonalMatrix([1,-1,-1,1]);
    ws_data[{1}][38] := DiagonalMatrix([1,-1,1,1]);

    return cover_data, ws_data;
end function;

procedure test_10_19()
    cover_data, ws_data := load_covers_and_ws_data_10_19();
    curves := GetHyperellipticCandidates();
    test_AllEquationsAboveCoversSingleCurve(10, 19, cover_data, ws_data, curves : manual_isomorphism);
    return;
end procedure;

test_10_19();