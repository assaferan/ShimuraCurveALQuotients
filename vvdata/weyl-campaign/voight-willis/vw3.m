// Voight-Willis (Computing power series expansions of modular forms, Sec. 3), with Magma's
// SingularValueDecomposition used correctly: it returns (Sigma, U, V) with the singular values on
// the DIAGONAL OF THE FIRST value in descending order and the right singular vectors as the ROWS
// of the third.  The earlier attempt read the diagonal of the second (an orthogonal matrix), so
// every singular value it reported was meaningless.
function DiscMatrix(mrep, q, p, CC)
    M := ChangeRing(mrep(q), CC);
    Phi := Matrix(CC, 2, 2, [1, -p, 1, -ComplexConjugate(p)]);
    W := Phi * M * Phi^(-1);
    return W / Sqrt(Determinant(W));
end function;
function DiscJ(W, w) return W[2][1]*w + W[2][2]; end function;

// real 2m x 2n form of a complex matrix: same singular values, each twice
function RealForm(Mx)
    CC := BaseRing(Mx); RR := RealField(Precision(CC));
    m := Nrows(Mx); n := Ncols(Mx);
    MR := Matrix(RR, 2*m, 2*n, [ RR | 0 : i in [1..4*m*n] ]);
    for i in [1..m] do for j in [1..n] do
        a := Real(Mx[i][j]); b := Imaginary(Mx[i][j]);
        MR[i][j] := a; MR[i][n+j] := -b; MR[m+i][j] := b; MR[m+i][n+j] := a;
    end for; end for;
    return MR;
end function;

// ascending singular values of a complex matrix (each appears twice; we de-duplicate by step 2)
function SingVals(Mx)
    MR := RealForm(Mx);
    Sig := SingularValueDecomposition(MR);
    r := Min(Nrows(MR), Ncols(MR));
    sv := Reverse([ Abs(Sig[i][i]) : i in [1..r] ]);      // ascending
    return [ sv[i] : i in [1..#sv by 2] ];                 // one per complex singular value
end function;

// coefficient vector of the smallest singular value
function SmallestVec(Mx)
    CC := BaseRing(Mx); n := Ncols(Mx);
    MR := RealForm(Mx);
    Sig, L, Rt := SingularValueDecomposition(MR);
    r := Min(Nrows(MR), Ncols(MR));
    j := r;                                                // descending order: last is smallest
    v := [ CC | Rt[j][i] + CC.1*Rt[j][n+i] : i in [1..n] ];
    return v, Abs(Sig[r][r]);
end function;

function ReducedSamples(G, O, DD, mrep, p, rad, Q, CC)
    out := []; gens := Generators(G);
    pi := Pi(RealField(Precision(CC)));
    for m in [1..Q] do
        wc := CC ! (rad * Exp(2*pi*CC.1*m/Q));
        red := ShimuraReduceUnit(G ! (O ! 1), gens, DD : z0 := DD ! wc);
        delta := red[1];
        if IsScalar(Quaternion(delta)) then continue; end if;
        W := DiscMatrix(mrep, Quaternion(delta), CC ! p, CC);
        wp := (W[1][1]*wc + W[1][2]) / (W[2][1]*wc + W[2][2]);
        Append(~out, <wc, wp, W>);
    end for;
    return out;
end function;

// eq (9) divided by (1-w')^k:  sum b_n (w'^n - jD^k w^n) = 0
function AutSystem(samples, k, nterms, CC)
    rows := [];
    for s in samples do
        w := s[1]; wp := s[2]; jk := DiscJ(s[3], w)^k;
        Append(~rows, [ wp^n - jk*w^n : n in [0..nterms-1] ]);
    end for;
    return Matrix(CC, #rows, nterms, &cat rows);
end function;

// f(z) = (1-w)^k sum b_n w^n, evaluated at a disc point w
function EvalForm(bs, k, w)
    return (1-w)^k * &+[ bs[n+1]*w^n : n in [0..#bs-1] ];
end function;

// Columns scale like rho^n, so without equilibration every singular value sits near rho^N and the
// gap is meaningless.  Scaling columns to unit norm corresponds to rescaling b_n.
function NormaliseColumns(Mx)
    for j in [1..Ncols(Mx)] do
        nm := Sqrt(&+[ Abs(Mx[i][j])^2 : i in [1..Nrows(Mx)] ]);
        if nm gt 0 then for i in [1..Nrows(Mx)] do Mx[i][j] := Mx[i][j]/nm; end for; end if;
    end for;
    return Mx;
end function;
