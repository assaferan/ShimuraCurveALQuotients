// The full m = 0-type sum over NONZERO cosets of the plane at X_0^15(2), with the fibre over each
// coset of L^v/L (prop:mult of the standalone, without the direct-sum hypothesis):
//
//     T = sum_{nu in L_-^v/L_-, nu != 0}  kappa^-_nu(0)  sum_{x in L_+^v : x + nu in L^v}  c_{[x+nu]}(-Q(x)),
//
// predicted correction per CM point = -T/4 (the fundamental case gives (1/2) c_eta(0) log N).
// kappa^-_nu(0) = log 2 * (W_nu/W_0)'(X = 1), with W_nu(X) = (1-X) sum_k alpha_k X^k the counting
// series of the coset at 2 (lem:pole), reconstructed as a rational function from alpha_0..alpha_K.
// c_eta(-m): principal part of F_f -- the cusp-0 coefficient of q^{-mM} at every eta with
// Q(eta) = m mod 1 (Guo-Yang Lemma 24), plus the cusp-oo coefficient of q^{-m} at eta = 0;
// c_eta(0) at a nonzero isotropic eta is 2 * M0MultiplierExact.
// Compare with fibreval.m (the code's value with every nonzero-coset m = 0-type term removed).
AttachSpec("ShimuraQuotients.spec");
SetColumns(0);
P<X> := PolynomialRing(Rationals());
D := 15; N := 2; KMAX := 10;
Ld := ShimuraCurveLattice(D, N);
Qg := ChangeRing(Ld`Q, Rationals());                 // Gram of <,> on L = Z^3; Q(x) = (x Qg, x)/2
n := Nrows(Qg);
e := func< i | Vector(Rationals(), [j eq i select 1 else 0 : j in [1..n]]) >;
Qf := func< v | (v*Qg, v)/2 >;
inLdual := func< v | forall{i : i in [1..n] | IsIntegral((v*Qg, e(i)))} >;
inL := func< v | forall{c : c in Eltseq(v) | IsIntegral(c)} >;
M := IsOdd(D*N) select 4*D*N else 2*D*N;

Xstar := CreateShimuraQuot(D, N, Set(Divisors(D*N)));
Xstar`g := GenusShimuraCurveQuotient(D, N, Xstar`W); Xstar`CurveID := 0;
curves := GetQuotientsAndGenera([Xstar]);
_ := exists(star){c : c in curves | IsStarCurve(c)};
fs := BorcherdsForms(star, curves : Prec := 100);
ks := Sort([k : k in Keys(fs)]);
etas := [fs[k] : k in ks];
mults := M0MultiplierExact(etas, Ld, D, N);
foos := [qExpansionAtoo(eta, 1) : eta in etas];
f0s := [qExpansionAt0(eta, 1) : eta in etas];
printf "forms %o, m = 0 multipliers (1/2)c_eta(0) = %o\n", ks, mults;
printf "pole orders at oo %o, at 0 (in q^(1/%o)) %o\n", [-Valuation(f) : f in foos], M, [-Valuation(f) : f in f0s];
bound := Maximum([-Valuation(f) : f in foos] cat [-Valuation(f)/M : f in f0s]);

// c_eta(-m) for the vector w in L^v representing eta, m = Q(x) > 0; c_eta(0) for eta nonzero isotropic
function coeff(i, w, m)
    if m eq 0 then
        assert not inL(w) and IsIntegral(Qf(w));
        return 2*mults[i];
    end if;
    c := Rationals()!0;
    if IsIntegral(m*M) then c +:= Coefficient(f0s[i], -Integers()!(m*M)); end if;
    if inL(w) and IsIntegral(m) then c +:= Coefficient(foos[i], -Integers()!m); end if;
    return c;
end function;

// rational reconstruction of sum alpha_k X^k from alpha_0..alpha_K: denominator degree r, numerator degree <= r
function pade(alphas, r)
    K := #alphas - 1;
    if K lt 2*r + 2 then return false, _; end if;
    // unknowns d_1..d_r; coefficient of X^t in D*A vanishes for t = r+1..K
    rows := [[alphas[t - j + 1] : j in [1..r]] : t in [r+1..K]];
    rhs := [-alphas[t + 1] : t in [r+1..K]];
    A := Matrix(Rationals(), rows); b := Vector(Rationals(), rhs);
    ok, sol := IsConsistent(Transpose(A), b);
    if not ok then return false, _; end if;
    Dn := 1 + &+[P | sol[j]*X^j : j in [1..r]];
    prod := Dn * &+[P | alphas[k+1]*X^k : k in [0..K]];
    Nm := &+[P | Coefficient(prod, k)*X^k : k in [0..r]];
    return true, Nm/Dn;
end function;

for d in [-15, -7, -52, -28, -60, -240, -48, -12] do
    lam := ChangeRing(ElementOfNorm(Ld`Q, -d, Ld`O, Ld`basis_L), Rationals());
    assert Qf(lam) eq -d;
    lam0 := lam / Content(ChangeRing(lam, Integers()));     // primitive vector on the CM line
    nn := Integers()!(lam0*Qg, lam0);                        // <lam0,lam0> = 2 Q(lam0)
    Mx := Matrix(Integers(), n, 1, [Integers() | (e(i)*Qg, lam) : i in [1..n]]);
    K := ChangeRing(KernelMatrix(Mx), Rationals());         // basis of L_-
    gram := K*Qg*Transpose(K);
    det := Integers()!Determinant(gram);
    idx := AbsoluteValue(Determinant(VerticalJoin(Matrix(1, n, Eltseq(lam0)), K)));
    d0 := FundamentalDiscriminant(d);
    printf "\n=== d = %o (d0 %o, conductor %o, 2 %o)  Q(lam0) = %o  [L : L_+ (+) L_-] = %o  det L_- = %o\n",
        d, d0, Isqrt(d div d0), KroneckerSymbol(d0, 2) eq 1 select "split" else "INERT", nn/2, idx, Factorization(det);
    // all cosets of L_-^v/L_- = Z^2 gram^-1 / Z^2: with P gram Q = S (Smith form) they are the
    // vectors (a/s1, b/s2) P, 0 <= a < s1, 0 <= b < s2
    dd := AbsoluteValue(det);
    Sm, Pm := SmithForm(ChangeRing(gram, Integers()));
    Pq := ChangeRing(Pm, Rationals());
    s1 := Sm[1,1]; s2 := Sm[2,2];
    cos := {};
    for a in [0..s1-1], b in [0..s2-1] do
        v := Vector(Rationals(), [a/s1, b/s2]) * Pq;
        key := [c - Floor(c) : c in Eltseq(v)];
        assert forall{c : c in Eltseq(v*gram) | IsIntegral(c)};
        Include(~cos, key);
    end for;
    assert #cos eq dd;
    integral := [key : key in cos | key ne [0,0] and IsIntegral((Vector(Rationals(), key)*gram, Vector(Rationals(), key))/2)];
    printf "  %o cosets, %o nonzero integral: %o\n", #cos, #integral, integral;
    // 2-adic counting series of a coset, A(X) = sum alpha_k X^k (alpha_0 = 1 for an integral coset)
    A11 := Integers()!(gram[1,1]/2); A12 := Integers()!gram[1,2]; A22 := Integers()!(gram[2,2]/2);
    series := function(key)
        a := Maximum([Valuation(Denominator(c), 2) : c in key] cat [0]);
        r1 := Integers()!(key[1]*2^a); r2 := Integers()!(key[2]*2^a);
        alphas := [Rationals()!1];
        for k in [1..KMAX] do
            cnt := 0; mod_ := 2^(k + 2*a); sh := 2^a;
            for u in [0..2^k-1] do
                U := r1 + sh*u;
                for t in [0..2^k-1] do
                    T := r2 + sh*t;
                    if (A11*U^2 + A12*U*T + A22*T^2) mod mod_ eq 0 then cnt +:= 1; end if;
                end for;
            end for;
            Append(~alphas, cnt / 2^k);
        end for;
        return alphas;
    end function;
    W := AssociativeArray();
    for key in [[0,0]] cat integral do
        if exists{c : c in key | Denominator(c) ne 2^Valuation(Denominator(c), 2)} then
            printf "  coset %o has odd part: skipped (anisotropic there)\n", key; continue;
        end if;
        alphas := series(key);
        ok := false;
        for r in [1..4] do
            ok, A := pade(alphas, r);
            if ok then break; end if;
        end for;
        error if not ok, Sprintf("no rational reconstruction for coset %o: %o", key, alphas);
        W[key] := (1 - X) * A;
        printf "  coset %-14o alphas %o  W = %o\n", key, alphas, W[key];
    end for;
    R0 := W[[0,0]];
    T := AssociativeArray(); for i in [1..#etas] do T[i] := Rationals()!0; end for;
    for key in integral do
        if not IsDefined(W, key) then continue; end if;
        R := W[key] / R0;
        error if Evaluate(R, 1) ne 0, Sprintf("ratio does not vanish at s = 0 for coset %o: %o", key, R);
        kap := Evaluate(Derivative(R), 1);                   // kappa^-_nu(0) / log 2
        nu := Vector(Rationals(), key) * K;
        for i in [1..#etas] do
            S := Rationals()!0; pairs := [];
            jmax := Floor(Sqrt(2*nn*bound)) + 1;
            for j in [-jmax..jmax] do
                x := j * lam0 / nn;                          // runs over L_+^v
                if not inLdual(x + nu) then continue; end if;
                m := Qf(x);
                if m gt bound then continue; end if;
                c := coeff(i, x + nu, m);
                S +:= c;
                if c ne 0 then Append(~pairs, <j, m, inL(x + nu) select "eta=0" else "eta!=0", c>); end if;
            end for;
            T[i] +:= kap * S;
            printf "  nu = %-14o form %-3o kappa/log2 = %-6o in L^v: %-5o sum_x c = %o  terms (j, Q(x), coset, c): %o\n",
                key, ks[i], kap, inLdual(nu), S, pairs;
        end for;
    end for;
    printf "  PREDICTED correction per point, -T/4, in units of log 2: %o   (fundamental recipe: %o)\n",
        [T[i]/-4 : i in [1..#etas]], mults;
end for;
