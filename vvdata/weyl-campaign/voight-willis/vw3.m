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

// A positive definite MAJORANT of the reduced norm, used wherever elements of an order (or of an
// ideal in it) of a prescribed reduced norm have to be enumerated.  The reduced norm form itself
// is indefinite of signature (2,2), so its level sets are hyperboloids and a cube in the
// coordinates of a basis is the wrong shape: it cuts them along directions unrelated to the
// lattice, and misses elements.  With iota the real matrix representation used throughout here,
// F(x) = Tr(iota(x) iota(x)^t) is positive definite, and for a real 2x2 matrix F >= 2|det|, so
// every element of reduced norm m has F(x) >= 2m and each class modulo the unit group has a
// shortest representative.  Enumerating F(x) <= C*m and raising C therefore terminates.
//
// The Gram matrix is real, so it is scaled by SCL and rounded; that can only change which
// candidates are OFFERED, never which are accepted, because every acceptance test downstream
// (reduced norm, membership in the order) is exact arithmetic.
SCL := 10^10;
function FrobeniusGram(bas, mrep, A)
    Ms := [ [ x : x in Eltseq(mrep(A ! b)) ] : b in bas ];
    n := #bas;
    GM := ZeroMatrix(Integers(), n, n);
    for i in [1..n] do for j in [i..n] do
        g := Round(SCL * (&+[ Ms[i][r]*Ms[j][r] : r in [1..4] ]));
        GM[i][j] := g; GM[j][i] := g;
    end for; end for;
    return GM;
end function;

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

// The side pairing of the fundamental domain, as SU(1,1) matrices with their inverses.  This is
// the set a reduction must descend over, and it is NOT Generators(G): that is a minimal generating
// set, three elements where the domain has fourteen sides at D = 15.  Magma records the side
// pairing in ShimFDGens, for the domain it most recently computed, so the domain for THIS disc is
// recomputed first.  Side pairings have reduced norm 1, so DiscMatrix divides by 1 and the matrix
// product below carries no square-root sign ambiguity.
function SidePairingMatrices(G, DD, mrep, p, CC)
    _ := FundamentalDomain(G, DD);
    try
        qs := [ Algebra(Parent(g)) ! g : g in G`ShimFDGens ];
    catch e
        qs := [ Quaternion(g) : g in Generators(G) ];
    end try;
    Ws := [];
    for q in qs do
        M := DiscMatrix(mrep, q, CC ! p, CC);
        Append(~Ws, M); Append(~Ws, M^(-1));
    end for;
    return Ws;
end function;

// Magma's fundamental domain here is the Dirichlet domain centred at the disc centre, so a point
// lies in it exactly when no side pairing decreases its modulus, and applying the pairing that
// decreases it most terminates inside.  Returns the accumulated Mobius matrix and the reduced
// point.
function ReduceByPairing(Ws, w, CC : maxit := 200)
    cur := w; bw := w;
    Acc := IdentityMatrix(CC, 2); bM := Acc;
    eps := (RealField(Precision(CC)) ! 10)^(-Precision(CC) div 3);
    for it in [1..maxit] do
        best := Abs(cur); moved := false;
        for M in Ws do
            jj := DiscJ(M, cur);
            if Abs(jj) le eps then continue; end if;
            nw := (M[1][1]*cur + M[1][2]) / jj;
            if Abs(nw) lt best - eps then best := Abs(nw); bM := M; bw := nw; moved := true; end if;
        end for;
        if not moved then break; end if;
        Acc := bM * Acc; cur := bw;
    end for;
    return Acc, cur;
end function;

// Sample a circle of radius rad and keep the points that carry a relation.  A point inside the
// domain reduces to itself and carries none, which is why such points are dropped; but
// ShimuraReduceUnit, given only Generators(G), returns the identity on points it FAILS to reduce
// as well, and those were being dropped with them.  At D = 15 that left 15% of a circle of radius
// rho, and 25% of a circle lying entirely outside the domain where the true figure must be 100%.
// The survivors cluster on short arcs, and fitting a degree-90 polynomial on a short arc is
// ill-conditioned however many points it has: that, not the size of the domain, is what made the
// larger bases fail.
function ReducedSamples(G, O, DD, mrep, p, rad, Q, CC)
    Ws := SidePairingMatrices(G, DD, mrep, p, CC);
    RR := RealField(Precision(CC));
    pi := Pi(RR);
    eps := (RR ! 10)^(-Precision(CC) div 3);
    out := [];
    for m in [1..Q] do
        wc := CC ! (rad * Exp(2*pi*CC.1*m/Q));
        W, wp := ReduceByPairing(Ws, wc, CC);
        if Abs(wp - wc) le eps then continue; end if;
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
