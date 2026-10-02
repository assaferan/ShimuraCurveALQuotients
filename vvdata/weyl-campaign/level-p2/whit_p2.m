// m = 0 local Whittaker factor of a p^2-scaled split binary lattice, by brute-force counting, to
// check the values Kudla-Yang Thm 4.3 gives with l_1 = l_2 = 2.
//   W_{0,p}(s, char(mu+L)) = (1-X) G(X),  alpha_k = p^(-k) #{x in (mu+L)/p^k L : Q(x) = 0 mod p^k},
//   Q(x,y) = p^2 (x^2 - y^2).  One representative per class; classes are (exact order, ord_p t_mu).
P<X> := PolynomialRing(Rationals());
KMAX := 4;
for p in [3, 5] do
    e1 := 1; e2 := -1;
    reps := AssociativeArray(); count := AssociativeArray();
    for a in [0..p^2-1], b in [0..p^2-1] do
        if (e1*a^2 + e2*b^2) mod p^2 ne 0 then continue; end if;
        ea := (a eq 0) select 0 else 2 - Valuation(a, p);
        eb := (b eq 0) select 0 else 2 - Valuation(b, p);
        tmu := (e1*a^2 + e2*b^2) div p^2;
        key := <Maximum(ea, eb), (tmu eq 0) select 99 else Valuation(tmu, p)>;
        if not IsDefined(reps, key) then reps[key] := <a, b>; count[key] := 0; end if;
        count[key] +:= 1;
    end for;
    printf "\np = %o : isotropic cosets of the p^2-scaled split plane\n", p;
    for key in Sort([k : k in Keys(reps)]) do
        a := reps[key][1]; b := reps[key][2];
        alphas := [];
        for k in [1..KMAX] do
            cnt := 0;
            for u in [0..p^k-1], v in [0..p^k-1] do
                if (e1*(a + p^2*u)^2 + e2*(b + p^2*v)^2) mod p^(k+2) eq 0 then cnt +:= 1; end if;
            end for;
            Append(~alphas, cnt / p^k);
        end for;
        W := (1 - X) * (1 + &+[P | alphas[k]*X^k : k in [1..KMAX]]);
        printf "  order p^-%o, ord_p t_mu = %-3o : %-4o cosets, alpha_0..%o = %o\n      (1-X)G = %o\n",
               key[1], key[2] eq 99 select "inf" else Sprint(key[2]), count[key], KMAX, [1] cat alphas, W;
    end for;
end for;
