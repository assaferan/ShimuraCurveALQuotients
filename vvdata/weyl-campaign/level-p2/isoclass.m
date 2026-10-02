// The finer prediction the level-p^2 multiplier rests on: the 3p^2-2p-1 nonzero isotropic cosets
// split into 2p(p-1) of exact order p^-2, 2(p-1) of order p^-1 with p | t_mu, and (p-1)^2 of order
// p^-1 with p not| t_mu (the last contributing nothing).  Classified here on the real lattices by
// the order of the coset in the discriminant group at p.
AttachSpec("ShimuraQuotients.spec");
for b in [[6,25],[6,49]] do
    D := b[1]; N := b[2]; p := [q : q in PrimeDivisors(N)][1];
    Ld := ShimuraCurveLattice(D, N);
    Q := ChangeRing(Ld`Q, Rationals()); dn := Ld`denom;
    cnt := AssociativeArray();
    for g in Ld`disc_grp do
        v := ChangeRing(g @@ Ld`to_disc, Rationals());
        r := (v*Q, v)/(2*dn^2);
        if r ne Floor(r) then continue; end if;          // isotropic only
        // the dual vector is v/dn, so its exact order at p is p^(v_p(dn) - min_i v_p(v_i))
        k := Maximum(0, Valuation(dn, p) - Minimum([Valuation(c, p) : c in Eltseq(v) | c ne 0] cat [Valuation(dn, p)]));
        if not IsDefined(cnt, k) then cnt[k] := 0; end if;
        cnt[k] +:= 1;
    end for;
    printf "\n%o_%o (p = %o): isotropic cosets by the p-part of their order in L^v/L\n", D, N, p;
    for k in Sort([k : k in Keys(cnt)]) do
        printf "   order p^%o : %o\n", k, cnt[k];
    end for;
    printf "   predicted: p^0 (zero coset and the anisotropic-part-free ones) 1, p^1 : %o, p^2 : %o\n",
           2*(p-1) + (p-1)^2, 2*p*(p-1);
end for;
