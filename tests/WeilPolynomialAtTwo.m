// Regression test for the p = 2 (2-rank) branch of IsHypWeilPolynomial, taken when g is not in
// [3..6] and 2 does not divide DN.
//
// The check: over F_2 a hyperelliptic curve of 2-rank f has f+1 ramification points, and
// J[2](F_2bar) is the degree-0 part of the permutation module on them. So the unit-root factor
// of the L-polynomial, reduced mod 2, is prod (t^d - 1)/(t - 1) over the Frobenius orbit sizes
// d. That factor is the first f+1 entries of Reverse(Coefficients(WeilPolynomial(X,2))). The
// remaining entries are even because they lie past the slope-0 segment.
//
// Before the fix the branch compared all 2g+1 coefficients against length-(f+1) table entries.
// It therefore rejected every curve that reached it, and it errored when f = 0 (no slope-0
// segment). The only recorded verdict through it was X_0(55,9)/<w9,w55> (CurveID 10146,
// "WeilPolynomial with p = 2"). That verdict was wrong.
//
// Why each expected value is trustworthy (each is known independently of this filter):
//   (55,9)/<9,55>  g=7  HYPERELLIPTIC: V3 descends (9 in W) and the quotient by V3 has genus 0,
//                  so 2g+2 = 16 fixed points (CheckModularNonALInvolutionModSym -> "V3 W1").
//   (95,1)/<1>     g=7  HYPERELLIPTIC: HyperellipticALInvolution (recorded, curve #11289).
//   X_0(97)        g=7  NOT hyperelliptic: Ogg's list of hyperelliptic X_0(N).
//   X_0(103)       g=8  NOT hyperelliptic: Ogg's list.
//   X_0(187)/w17   g=7  2-rank 0: the f = 0 case, which used to error.
// Measured on 2026-09-25 over data/curves_after_UpdateCurves8.dat with 2 not dividing DN: the
// fixed p = 2 check passed all 337 decided-hyperelliptic curves (g = 2..9) and rejected 297 of
// 1120 decided non-hyperelliptic ones (g = 3..11).

function wpat2_setup(g)
    at2 := AssociativeArray();
    for f in [0..g] do at2[f] := HyperellipticWeilPolysAtTwo(f); end for;
    pw := AssociativeArray(); pw[g] := AssociativeArray();   // no odd primes: bound 2 isolates p = 2
    return pw, at2;
end function;

// id is the CurveID in data/curves_after_*.dat. It must be set: WeilPolynomial's point-count
// cache compares curves with 'eq', which reads CurveID.
function wpat2_curve(D, N, W, g, id)
    X := CreateShimuraQuot(D, N, W);
    X`g := g;
    X`CurveID := id;
    return X;
end function;

function wpat2_rank(X)
    slopes := SlopesWithMultiplicities(NewtonPolygon(WeilPolynomial(X, 2), 2));
    return &+[Integers() | s[2] : s in slopes | s[1] eq 0];
end function;

procedure test_WeilPolynomialAtTwo()
    assert HyperellipticWeilPolysAtTwo(0) eq [[GF(2) | 1]];

    // Hyperelliptic: must pass.
    for data in [<55, 9, {1, 9, 55, 495}, 7, 7, 10146>, <95, 1, {1}, 7, 7, 11287>] do
        D, N, W, g, f, id := Explode(data);
        X := wpat2_curve(D, N, W, g, id);
        assert wpat2_rank(X) eq f;
        pw, at2 := wpat2_setup(g);
        assert IsHypWeilPolynomial(X, pw, at2, 2);
    end for;

    // Non-hyperelliptic, correctly rejected at p = 2.
    for data in [<1, 97, {1}, 7, 7, 451>, <1, 103, {1}, 8, 8, 482>] do
        D, N, W, g, f, id := Explode(data);
        X := wpat2_curve(D, N, W, g, id);
        assert wpat2_rank(X) eq f;
        pw, at2 := wpat2_setup(g);
        b, p := IsHypWeilPolynomial(X, pw, at2, 2);
        assert not b and p eq 2;
    end for;

    // 2-rank 0: the unit-root factor is 1, so nothing is constrained and it must not error.
    X := wpat2_curve(1, 187, {1, 17}, 7, 868);
    assert wpat2_rank(X) eq 0;
    pw, at2 := wpat2_setup(7);
    assert IsHypWeilPolynomial(X, pw, at2, 2);
end procedure;

test_WeilPolynomialAtTwo();
