// Normalisation of the repo's local Whittaker polynomials at a CONDUCTOR prime, against a density
// count.  For a coset mu of L_- and m > 0 put alpha_k = p^-k #{x in (mu+L_-)/p^k L_- : Q(x) = m mod p^k}
// and W^count(X) = (1-X) sum_k alpha_k X^k; at a unimodular prime with p not dividing m this is
// 1 - chi(p) X/p, the Euler factor of 1/L(s+1,chi).  The repo's LocalWhittakerPolynomial (bare) times
// p^(-v_p(det S)/2) is Kudla-Yang's W without the Weil index.  Compare the two on X_0^15(2).
AttachSpec("ShimuraQuotients.spec");
SetColumns(0);
P<X> := PolynomialRing(Rationals());
D := 15; N := 2; KMAX := 9;
Ld := ShimuraCurveLattice(D, N);
Q := ChangeRing(Ld`Q, Integers()); Qr := ChangeRing(Q, Rationals());
function count_series(gram, key, m, p)
    KMAX := p eq 2 select 9 else 4;          // p^(2k) iterations per term
    A11 := gram[1,1]/2; A12 := gram[1,2]; A22 := gram[2,2]/2;
    a := Maximum([Valuation(Denominator(c), p) : c in key] cat [0]);
    r1 := Integers()!(key[1]*p^a); r2 := Integers()!(key[2]*p^a);
    // alpha_0 = 1 iff Q(mu) = m mod Z, i.e. the coset can represent m at all
    qmu := (Vector(Rationals(), key)*gram, Vector(Rationals(), key))/2;
    alphas := [Rationals() | IsIntegral(qmu - m) select 1 else 0];
    for k in [1..KMAX] do
        cnt := 0; md := p^(k + 2*a); sh := p^a; target := Integers()!(m*p^(2*a));
        for u in [0..p^k-1] do U := r1 + sh*u;
            for t in [0..p^k-1] do T := r2 + sh*t;
                if (Integers()!(A11*U^2 + A12*U*T + A22*T^2) - target) mod md eq 0 then cnt +:= 1; end if;
            end for;
        end for;
        Append(~alphas, cnt / p^k);
    end for;
    return (1 - X) * &+[P | alphas[k+1]*X^k : k in [0..KMAX]];
end function;
for d in [-15, -60, -240, -28, -12, -48] do
    lam := ElementOfNorm(Q, -d, Ld`O, Ld`basis_L);
    lam := ChangeRing(lam, Integers());
    Lminus := Kernel(Transpose(Matrix(lam*Q)));
    B := BasisMatrix(Lminus); S := B*Q*Transpose(B); gram := ChangeRing(S, Rationals());
    d0 := FundamentalDiscriminant(d); _, f := IsSquare(d div d0);
    printf "\n=== d = %o (d0 %o, f %o) det S = %o\n", d, d0, f, Factorization(Integers()!Determinant(S));
    for p in [2, 7] do
        v := Valuation(Integers()!Determinant(S), p);
        chi := KroneckerSymbol(d0, p);
        for m in [1, 2, 3, 4, 5, 6] do
            mu := Vector(Rationals(), [0, 0]) * ChangeRing(B, Rationals());   // the zero coset, as a vector in L
            bare := LocalWhittakerPolynomial(Rationals()!m, p, mu, Lminus, Q);
            cnt := count_series(gram, [0, 0], m, p);
            KM := p eq 2 select 9 else 4;
            trunc := func< F | &+[P | Coefficient(F, k)*X^k : k in [0..KM-1]] >;
            ratio := (trunc(cnt) eq 0 or trunc(P!bare) eq 0) select "n/a" else Sprint(trunc(cnt) / trunc(P!bare));
            printf "  p = %o (v_p det = %o, chi(p) = %2o) m = %o: repo bare %-30o count %-34o count/bare = %o\n",
                p, v, chi, m, bare, trunc(cnt), ratio;
        end for;
    end for;
end for;
