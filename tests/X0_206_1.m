import "tests/BorcherdsProducts.m" : test_AllEquationsAboveCoversSingleCurve;

function load_covers_and_ws_data_206_1()
    _<s> := PolynomialRing(Rationals());

    // D = 206
    cover_data := AssociativeArray();
    cover_data[{1}] := <HyperellipticCurve(-8*s^20+13*s^18+42*s^16+331*s^14+220*s^12-733*s^10-6646*s^8-19883*s^6-28840*s^4-18224*s^2-4096), Matrix([[0,0,-1],[0,64,0],[-1,0,0]])>;

    // ⚠ THE OTHER THREE COVER KEYS, added 2026-09-09. This test used to carry ONLY W={1}, so it
    // checked 1 of the model's 4 covers and silently skipped the rest -- the helper's
    // `if not is_def then continue` makes an absent key invisible rather than an error.
    // ⚠ [1,103] IS STORED AS <genus, f, h>, i.e. y^2 + h*y = f with h = x^5+x^4+x^3+x^2. Passing
    // only f would compare a DIFFERENT curve of the same genus -- exactly the defect that made
    // tests/_offline/X0_87_1.m fail for two days.
    cover_data[{1,2}] := <HyperellipticCurve(Polynomial(Rationals(), [ -1/512, 13/4096, 21/2048, 331/4096, 55/1024, -733/4096, -3323/2048, -19883/4096, -3605/512, -1139/256, -1 ])), DiagonalMatrix([1,1,1])>;   // genus 4
    cover_data[{1,206}] := <HyperellipticCurve(Polynomial(Rationals(), [ 0, 1 ])), DiagonalMatrix([1,1,1])>;   // genus 0
    cover_data[{1,103}] := <HyperellipticCurve(Polynomial(Rationals(), [ 0, 1024, -4556, 7210, -4971, 1661, -184, -56, 82, -11, 3, 2 ]), Polynomial(Rationals(), [ 0, 0, 1, 1, 1, 1 ])), DiagonalMatrix([1,1,1])>;   // genus 5, y^2 + h*y = f

    ws_data := AssociativeArray();
    ws_data[{1}] := AssociativeArray();
    ws_data[{1}][2] := DiagonalMatrix([-1,1,1]);
    ws_data[{1}][206] := DiagonalMatrix([1,-1,1]);

    return cover_data, ws_data;
end function;

procedure test_206_1()
    cover_data, ws_data := load_covers_and_ws_data_206_1();
    curves := GetHyperellipticCandidates();
    test_AllEquationsAboveCoversSingleCurve(206, 1, cover_data, ws_data, curves);
    return;
end procedure;

test_206_1();