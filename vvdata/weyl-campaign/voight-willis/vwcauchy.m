// The Cauchy-integral formulation of Voight-Willis, Sec. 3 (eq. 13-14) with the two rescalings of
// Remark 3.4.  This is what the paper's own Shimura-curve examples use; the raw automorphy system
// of eq. (9)-(11), which vw3.m implements, is presented there first but then superseded.
//
// Writing F(w) = sum b_n w^n for the bare series, the validated automorphy relation of AutSystem is
// F(w') = J F(w) with J = DiscJ(delta, w)^k, so F(w) = J^-1 F(w').  Cauchy on the circle |w| = R,
// discretised at Q equally spaced points, gives
//
//     b_n  =  (1/Q) sum_m F(w_m) / w_m^n  =  (1/Q) sum_m J_m^-1 sum_r b_r (w'_m)^r / w_m^n
//
// so b is in the kernel of K - 1 with K_{nr} = (1/Q) sum_m J_m^-1 (w'_m)^r / w_m^n.  The matrix is
// SQUARE, (N+1) x (N+1), where the automorphy system is overdetermined; and the trapezoidal rule on
// a periodic analytic integrand aliases as R^Q, so Q of the order of N already puts the quadrature
// error at the truncation floor.
//
// Remark 3.4's two rescalings are what make it well conditioned, and both matter:
//   * b'_n = b_n R^n, which replaces w_m by w_m/R of absolute value 1 and w'_m by w'_m/R of
//     absolute value at most 1.  Without it the entries grow like R^-n.
//   * divide row n by (K_nn - 1), so the diagonal of K - 1 is 1.  The paper observes that the
//     largest entry of each row is then the diagonal one, i.e. the matrix is diagonally dominant.

// All Q points of the circle of radius R, each with the element that reduces it into the domain
// (the identity when it is already inside).  Unlike ReducedSamples nothing is skipped: the Cauchy
// sum is over the whole circle.
function CauchySamples(G, O, DD, mrep, p, R, Q, CC)
    out := []; gens := Generators(G);
    pi := Pi(RealField(Precision(CC)));
    for m in [1..Q] do
        wc := CC ! (R * Exp(2*pi*CC.1*m/Q));
        red := ShimuraReduceUnit(G ! (O ! 1), gens, DD : z0 := DD ! wc);
        delta := red[1];
        W := DiscMatrix(mrep, Quaternion(delta), CC ! p, CC);
        wp := (W[1][1]*wc + W[1][2]) / (W[2][1]*wc + W[2][2]);
        Append(~out, <wc, wp, W>);
    end for;
    return out;
end function;

// K - 1, rescaled by Remark 3.4.  nterms = N+1 coefficients b_0..b_N.
function CauchySystem(samples, k, nterms, R, CC)
    Q := #samples;
    // powers, already divided by R: |what| <= 1 for both
    wh  := [ s[1]/R : s in samples ];           // modulus exactly 1
    whp := [ s[2]/R : s in samples ];           // modulus at most rho/R
    Jin := [ DiscJ(s[3], s[1])^(-k) : s in samples ];
    // pw[m][r+1] = (w'_m/R)^r ;  iw[m][n+1] = (w_m/R)^-n
    pw := [ [ CC | 1 : r in [1..nterms] ] : m in [1..Q] ];
    iw := [ [ CC | 1 : n in [1..nterms] ] : m in [1..Q] ];
    for m in [1..Q] do
        for r in [2..nterms] do pw[m][r] := pw[m][r-1]*whp[m]; end for;
        for n in [2..nterms] do iw[m][n] := iw[m][n-1]/wh[m]; end for;
    end for;
    A := ZeroMatrix(CC, nterms, nterms);
    for n in [1..nterms] do
        for r in [1..nterms] do
            A[n][r] := (&+[ Jin[m]*pw[m][r]*iw[m][n] : m in [1..Q] ]) / Q;
        end for;
        A[n][n] := A[n][n] - 1;                  // K - 1
    end for;
    for n in [1..nterms] do                      // row scaling: unit diagonal
        d := A[n][n];
        if Abs(d) gt 0 then for r in [1..nterms] do A[n][r] := A[n][r]/d; end for; end if;
    end for;
    return A;
end function;
