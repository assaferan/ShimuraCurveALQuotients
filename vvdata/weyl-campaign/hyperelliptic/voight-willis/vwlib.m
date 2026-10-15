// Voight-Willis style expansion layer (reconstructed from the method's description, not from the
// paper itself).  Expand a weight-k form in the unit-disc coordinate centred at the fundamental
// domain's centre, impose f(gw) = j(g,w)^k f(w) for each side-pairing generator g at sample points
// on g's isometric circle (where |j| = 1, so a point and its image both stay inside the domain
// radius and the truncated series is accurate at both), and read the space of forms off the
// numerical nullspace.

function DiscMatrix(mrep, q, p, CC)
    M := ChangeRing(mrep(q), CC);
    Phi := Matrix(CC, 2, 2, [1, -p, 1, -ComplexConjugate(p)]);
    W := Phi * M * Phi^(-1);
    return W / Sqrt(Determinant(W));
end function;

function DiscAct(W, w) return (W[1][1]*w + W[1][2]) / (W[2][1]*w + W[2][2]); end function;
function DiscJ(W, w)  return W[2][1]*w + W[2][2]; end function;

// Points on the isometric circle |cw+d| = 1 of W: centre -d/c, radius 1/|c|.
function IsoCirclePoints(W, npts, CC)
    c := W[2][1]; d := W[2][2];
    if Abs(c) lt 10^(-Precision(CC) div 2) then return []; end if;
    ctr := -d/c; rad := 1/Abs(c);
    pi := Pi(RealField(Precision(CC)));
    return [ ctr + rad*Exp(2*pi*CC.1*j/npts) : j in [0..npts-1] ];
end function;

// Linear system of the automorphy relations, truncating f = sum_{n<nterms} a_n w^n.
// Sample only the arc of each isometric circle that is an actual side of the fundamental domain,
// found by keeping the points whose image also stays within the domain radius: there both w and
// g w have modulus at most R, so the truncation error is of order R^nterms at every entry.
function AutomorphySystem(Ws, k, nterms, npts, R, CC)
    rows := []; maxmod := 0; used := 0;
    for W in Ws do
        for w in IsoCirclePoints(W, npts, CC) do
            if Abs(w) gt R then continue; end if;
            gw := DiscAct(W, w);
            if Abs(gw) gt R then continue; end if;
            used +:= 1;
            maxmod := Max([maxmod, Abs(w), Abs(gw)]);
            jk := DiscJ(W, w)^k;
            Append(~rows, [ gw^n - jk*w^n : n in [0..nterms-1] ]);
        end for;
    end for;
    if #rows eq 0 then return ZeroMatrix(CC, 1, nterms), 0.0, 0; end if;
    return Matrix(CC, #rows, nterms, &cat rows), maxmod, used;
end function;

// Scale every column to unit norm.  Without this the trailing columns, whose entries are of
// order R^n, are numerically negligible and contribute spurious near-null directions whose number
// grows with the truncation length rather than reflecting the space of forms.
function NormaliseColumns(Mx)
    CC := BaseRing(Mx);
    for j in [1..Ncols(Mx)] do
        nm := Sqrt(&+[ Abs(Mx[i][j])^2 : i in [1..Nrows(Mx)] ]);
        if nm gt 0 then for i in [1..Nrows(Mx)] do Mx[i][j] := Mx[i][j]/nm; end for; end if;
    end for;
    return Mx;
end function;

// Singular values, ascending.  The nullity is the count below a gap in this list.
function SingVals(Mx)
    CC := BaseRing(Mx);
    R := RealField(Precision(CC));
    m := Nrows(Mx); n := Ncols(Mx);
    // real 2m x 2n form of the complex matrix: same singular values, each doubled
    MR := Matrix(R, 2*m, 2*n, [ R | 0 : i in [1..4*m*n] ]);
    for i in [1..m] do for j in [1..n] do
        a := Real(Mx[i][j]); b := Imaginary(Mx[i][j]);
        MR[i][j] := a;       MR[i][n+j] := -b;
        MR[m+i][j] := b;     MR[m+i][n+j] := a;
    end for; end for;
    _, S, _ := SingularValueDecomposition(MR);
    sv := Sort([ Abs(S[i][i]) : i in [1..Min(Nrows(S), Ncols(S))] ]);
    return sv;
end function;
