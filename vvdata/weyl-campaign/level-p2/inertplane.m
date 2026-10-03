// The m = 0 local factor at 2 on the ACTUAL negative plane L_- = L cap lambda^perp of X_0^15(2),
// for the conductor discriminants, by brute-force counting (the normalisation of lem:m0loc):
//   W_{0,2}(s, char(mu+L)) = (1-X) G(X),  alpha_k = 2^(-k) #{x in (mu+L_-)/2^k L_- : Q(x) = 0 mod 2^k}.
// Question: does the ZERO coset have a simple pole (alpha_k growing linearly) at the 4-scaled
// inert plane (d = -48), as it does at every split plane, or no pole, as at the 2-scaled inert
// plane (d = -12)?  The whole m = 0 correction at the level prime rests on that pole.
AttachSpec("ShimuraQuotients.spec");
P<X> := PolynomialRing(Rationals());
D := 15; N := 2; KMAX := 7;
Ld := ShimuraCurveLattice(D, N);
Q := ChangeRing(Ld`Q, Rationals());
e := func< i | Vector(Rationals(), [j eq i select 1 else 0 : j in [1..3]]) >;
for d in [-15, -60, -240, -12, -48] do
    lam := ElementOfNorm(Ld`Q, -d, Ld`O, Ld`basis_L);
    lv := ChangeRing(lam, Rationals());
    Mx := Matrix(Integers(), 3, 1, [Integers() | (e(i)*Q, lv) : i in [1..3]]);
    K := KernelMatrix(Mx);                     // basis of L_- (rank 2)
    gram := Matrix(Rationals(), 2, 2, [(ChangeRing(K[i], Rationals())*Q, ChangeRing(K[j], Rationals())) : i, j in [1..2]]);
    G := ChangeRing(gram, Integers());
    det := Determinant(G);
    d0 := FundamentalDiscriminant(d);
    printf "\n=== d = %o (fundamental %o, conductor %o, 2 %o in K)  L_- Gram = %o  det %o = 2^%o * %o\n",
        d, d0, Integers()!Sqrt(d/d0), KroneckerSymbol(d0, 2) eq 1 select "split" else (KroneckerSymbol(d0,2) eq -1 select "INERT" else "ramified"),
        Eltseq(G), det, Valuation(det, 2), det div 2^Valuation(det, 2);
    // cosets of the 2-part of L_-^v / L_-: v = w / 2^a with 2^a || det, v*G integral
    a := Valuation(det, 2); dd := 2^a;
    seen := {};
    for w in CartesianPower([0..dd-1], 2) do
        v := Vector(Rationals(), [w[1]/dd, w[2]/dd]);
        if not forall{c : c in Eltseq(v*gram) | IsIntegral(c)} then continue; end if;
        key := [c - Floor(c) : c in Eltseq(v)];
        if key in seen then continue; end if;
        Include(~seen, key);
        r := (v*gram, v)/2;
        alphas := [];
        for k in [1..KMAX] do
            cnt := 0;
            for u, t in [0..2^k-1] do
                x := Vector(Rationals(), [key[1] + u, key[2] + t]);
                q := (x*gram, x)/2;
                if IsIntegral(q) and (Integers()!q) mod 2^k eq 0 then cnt +:= 1; end if;
            end for;
            Append(~alphas, cnt / 2^k);
        end for;
        W := (1 - X) * (1 + &+[P | alphas[k]*X^k : k in [1..KMAX]]);
        pole := (alphas[KMAX] gt alphas[KMAX-2]) select "grows: POLE" else "bounded: no pole";
        printf "  mu = %-14o Q(mu) = %-6o %-12o alphas %o  %o\n", key, r,
            (r eq Floor(r)) select "integral" else "non-integral", alphas, pole;
    end for;
end for;
