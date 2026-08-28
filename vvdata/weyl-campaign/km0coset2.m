// Same test, restricted to the cosets that actually occur: mu in L^dual/L.
// For A = p*(unimodular), L^dual/L = (Z/p)^2 and the nonzero cosets are
// mu = (a/p, b/p).  For A unimodular, L^dual/L = 0 -- nothing to test, which is
// exactly the regime Schofer's cited Thm 5.2 was proved in.
AttachSpec("ShimuraQuotients.spec");
L := RSpaceWithBasis(IdentityMatrix(Integers(), 2));
z := Rationals()!0;
for p in [2,3,5,7] do
    u := 1; while KroneckerSymbol(u,p) ne -1 do u +:= 1; end while;
    for t -> A in [* Matrix(Integers(),2,2,[0,p,p,0]),
                     Matrix(Integers(),2,2,[2*p,0,0,2*p*u]) *] do
        lab := t eq 1 select "p-scaled hyperbolic (p SPLIT)" else "p-scaled diagonal  (p inert)";
        nz := [];
        for a in [0..p-1] do for b in [0..p-1] do
            if a eq 0 and b eq 0 then continue; end if;
            w := LocalWhittakerAtOne(z, p, Vector([Rationals()|a/p, b/p]), L, A);
            if w ne 0 then Append(~nz, <a, b, w>); end if;
        end for; end for;
        printf "p=%-2o %o : W_{0,p} != 0 at %o of %o nonzero cosets of L^dual/L%o\n",
            p, lab, #nz, p^2-1,
            #nz eq 0 select "" else Sprintf("   values %o", {* x[3] : x in nz *});
    end for;
end for;
printf "DONE\n"; quit;
