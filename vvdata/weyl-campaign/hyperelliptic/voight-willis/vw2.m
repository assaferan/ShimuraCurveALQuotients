// Voight-Willis as published (Computing power series expansions of modular forms, Sec. 3).
// f(z) = (1-w)^k sum b_n w^n;  sample w on the circle |w| = rho, REDUCE it into the Dirichlet
// domain with g so that w' = gw in D, and impose eq (9).  Dividing eq (9) by (1-w')^k the relation
// is sum b_n (w'^n - jD^k w^n) = 0 with jD the disc automorphy factor (checked equal to
// j(g,z)^k (1-w)^k/(1-w')^k to 1e-60).  The earlier attempt used the same relation but sampled on
// the generators' isometric circles instead; that is the only change.
function DiscMatrix(mrep, q, p, CC)
    M := ChangeRing(mrep(q), CC);
    Phi := Matrix(CC, 2, 2, [1, -p, 1, -ComplexConjugate(p)]);
    W := Phi * M * Phi^(-1);
    return W / Sqrt(Determinant(W));
end function;
function DiscJ(W, w) return W[2][1]*w + W[2][2]; end function;

function NormaliseColumns(Mx)
    for j in [1..Ncols(Mx)] do
        nm := Sqrt(&+[ Abs(Mx[i][j])^2 : i in [1..Nrows(Mx)] ]);
        if nm gt 0 then for i in [1..Nrows(Mx)] do Mx[i][j] := Mx[i][j]/nm; end for; end if;
    end for;
    return Mx;
end function;

function SingVals(Mx)
    CC := BaseRing(Mx); RR := RealField(Precision(CC));
    m := Nrows(Mx); n := Ncols(Mx);
    MR := Matrix(RR, 2*m, 2*n, [ RR | 0 : i in [1..4*m*n] ]);
    for i in [1..m] do for j in [1..n] do
        a := Real(Mx[i][j]); b := Imaginary(Mx[i][j]);
        MR[i][j] := a; MR[i][n+j] := -b; MR[m+i][j] := b; MR[m+i][n+j] := a;
    end for; end for;
    _, S, _ := SingularValueDecomposition(MR);
    return Sort([ Abs(S[i][i]) : i in [1..Min(Nrows(S), Ncols(S))] ]);
end function;

// Sample points of modulus rho, each reduced into the domain; returns pairs <w, w', W_disc>.
function ReducedSamples(G, O, DD, mrep, p, rho, Q, CC)
    out := []; gens := Generators(G);
    pi := Pi(RealField(Precision(CC)));
    for m in [1..Q] do
        wc := CC ! (rho * Exp(2*pi*CC.1*m/Q));
        wpt := DD ! wc;
        red := ShimuraReduceUnit(G ! (O ! 1), gens, DD : z0 := wpt);
        delta := red[1];
        if IsScalar(Quaternion(delta)) then continue; end if;
        W := DiscMatrix(mrep, Quaternion(delta), CC ! p, CC);
        wp := (W[1][1]*wc + W[1][2]) / (W[2][1]*wc + W[2][2]);
        Append(~out, <wc, wp, W>);
    end for;
    return out;
end function;

// Paper Remark 3.3: work with f_norm(z) = f(z) (Im z)^{k/2}, whose automorphy factor
// (j/|j|)^k has absolute value 1.  With f = (1-w)^k sum b_n w^n the relation becomes
//   (Im gz)^{k/2} (1-w')^k sum b_n w'^n = (j/|j|)^k (Im z)^{k/2} (1-w)^k sum b_n w^n,
// every factor of size O(1).  Without this the rows carry |j|^k, which swamps the residual.
function AutSystemNorm(samples, k, nterms, p, CC)
    rows := [];
    pc := CC ! p;
    zOf := func< w | (ComplexConjugate(pc)*w - pc)/(w - 1) >;
    for s in samples do
        w := s[1]; wp := s[2]; jD := DiscJ(s[3], w);
        u := (jD/Abs(jD))^k;
        yz  := Imaginary(zOf(w));
        ygz := Imaginary(zOf(wp));
        ca := ygz^(k div 2) * (1-wp)^k;
        cb := u * yz^(k div 2) * (1-w)^k;
        Append(~rows, [ ca*wp^n - cb*w^n : n in [0..nterms-1] ]);
    end for;
    return Matrix(CC, #rows, nterms, &cat rows);
end function;

function RowNormalise(Mx)
    for i in [1..Nrows(Mx)] do
        nm := Sqrt(&+[ Abs(Mx[i][j])^2 : j in [1..Ncols(Mx)] ]);
        if nm gt 0 then for j in [1..Ncols(Mx)] do Mx[i][j] := Mx[i][j]/nm; end for; end if;
    end for;
    return Mx;
end function;

// Plain disc form of eq (9) after dividing by (1-w')^k: sum b_n (w'^n - jD^k w^n) = 0.
function AutSystem(samples, k, nterms, CC)
    rows := [];
    for s in samples do
        w := s[1]; wp := s[2]; jk := DiscJ(s[3], w)^k;
        Append(~rows, [ wp^n - jk*w^n : n in [0..nterms-1] ]);
    end for;
    return Matrix(CC, #rows, nterms, &cat rows);
end function;
