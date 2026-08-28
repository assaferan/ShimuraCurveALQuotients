// THE DECISIVE LOCAL TEST.  Schofer Lemma 2.18 kills the m=0 term at every
// coset phi != phi_0, citing "Thm 5.2 of [16]: W_{0,p}(s,phi_p) = 0 if phi_p is
// not the characteristic function of the local lattice."  Does that hold when
// the local lattice is NOT unimodular -- i.e. at a split level prime where L_-
// is N-SCALED?  If W_{0,p} is nonzero at some mu != 0 there, that IS the dropped
// term, and Lemma 2.18's proof is being applied outside its hypotheses.
AttachSpec("ShimuraQuotients.spec");
L := RSpaceWithBasis(IdentityMatrix(Integers(), 2));
zero := Rationals()!0;
for p in [2,3,5] do
    u := 1; while KroneckerSymbol(u,p) ne -1 do u +:= 1; end while;
    for sc in [1, p] do
      for nm -> A in [* Matrix(Integers(),2,2,[0,sc,sc,0]),
                        Matrix(Integers(),2,2,[2*sc,0,0,2*sc*u]) *] do
        lab := nm eq 1 select "hyperbolic" else "diagonal  ";
        // cosets mu in (1/(p*sc)) L / L  -- enough to see the scaled dual
        nz := []; 
        den := p*sc;
        for a in [0..den-1] do for b in [0..den-1] do
            if a eq 0 and b eq 0 then continue; end if;
            mu := Vector([Rationals() | a/den, b/den]);
            w := LocalWhittakerAtOne(zero, p, mu, L, A);
            if w ne 0 then Append(~nz, <a, b, w>); end if;
        end for; end for;
        printf "p=%o scale=%-2o %o : W_{0,p} nonzero at %o of %o nonzero cosets%o\n",
            p, sc, lab, #nz, den^2-1,
            #nz eq 0 select "" else Sprintf("   e.g. %o", nz[1..Minimum(3,#nz)]);
      end for;
    end for;
end for;
printf "DONE\n"; quit;
