// km0coset2.m measured the VALUE W_{0,p}(0,mu) at the nonzero cosets and found it
// equal to 1 on exactly the 2(p-1) isotropic ones.  That is not enough to decide the
// provenance of the log N.
//
// Schofer Thesis, Lemma 2.18 proof, verbatim:
//     E_0(tau,s;phi,+1) = v^(s/2) phi(0)
//                         - sqrt(pi) i v^(-s/2) Gamma((s+1)/2)/Gamma(s/2+1)
//                           * prod_p W_{0,p}(s,phi_p).
// For mu != 0 we have phi_mu(0) = 0, so the first term drops; the archimedean factor
// at s=0 is -i*pi, NOT zero; and the value vanishes by incoherence.  Hence
//     kappa_mu(0) = b_{phi_mu}(0) = -i*pi * d/ds [ prod_p W_{0,p}(s,phi_mu) ]_{s=0}.
// A log N can therefore enter only if the LEVEL-PRIME factor is a nonconstant
// polynomial in X = N^(-s).  Its value being 1 is compatible with either.
//
// So: is W_{0,N}(X,mu) identically 1, or nonconstant with W(1) = 1?
AttachSpec("ShimuraQuotients.spec");
L := RSpaceWithBasis(IdentityMatrix(Integers(), 2));
z := Rationals()!0;
for p in [2,3,5,7] do
    u := 1; while KroneckerSymbol(u,p) ne -1 do u +:= 1; end while;
    for t -> A in [* Matrix(Integers(),2,2,[0,p,p,0]),
                     Matrix(Integers(),2,2,[2*p,0,0,2*p*u]) *] do
        lab := t eq 1 select "p-scaled hyperbolic (p SPLIT)" else "p-scaled diagonal  (p inert)";
        printf "\n===== p = %o : %o =====\n", p, lab;
        polys := {* *};
        nz := 0;
        for a in [0..p-1] do for b in [0..p-1] do
            if a eq 0 and b eq 0 then continue; end if;
            mu := Vector([Rationals()| a/p, b/p]);
            W := LocalWhittakerPolynomial(z, p, mu, L, A);
            v1 := Evaluate(W, 1);
            if v1 ne 0 then
                nz +:= 1;
                Include(~polys, W);
                if nz le 3 then
                    printf "   mu=(%o/%o,%o/%o)  W(X) = %-24o W(1) = %-6o W'(1) = %o\n",
                           a, p, b, p, W, v1, Evaluate(Derivative(W), 1);
                end if;
            end if;
        end for; end for;
        printf "   %o of %o nonzero cosets have W(1) != 0   (2(p-1) = %o)\n", nz, p^2-1, 2*(p-1);
        printf "   DISTINCT polynomials among them: %o\n", {@ f : f in polys @};
        allconst := &and[ Degree(f) le 0 : f in polys ];
        printf "   >>> level-prime factor is %o  ==>  %o\n",
               allconst select "CONSTANT in X" else "NONCONSTANT in X",
               allconst select "NO log p can come from this place"
                            else "a log p CAN come from this place";
    end for;
end for;
printf "\nDONE\n"; quit;
