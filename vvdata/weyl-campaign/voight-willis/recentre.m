// Magma's UnitDisc defaults to Center := (9/10)*i, a fixed point of the upper half plane with no
// relation to the group.  The expansion is centred there, so the circumradius rho of the
// fundamental domain -- which sets both the number of terms (rho^N) and the conditioning -- is an
// accident of that choice.  Recentring on the domain's own centre shrinks rho.
//
// Iterate: compute the domain, take the Euclidean centroid of its vertices in the disc, map that
// back to the upper half plane and rebuild the disc there.  Report rho at each step.
//
//   magma -b DISC:=10 LEVEL:=1 PR:=100 ITERS:=4 recentre.m < /dev/null
AttachSpec("ShimuraQuotients.spec");
SetColumns(0);
D := StringToInteger(DISC); N := StringToInteger(LEVEL);
PREC := StringToInteger(PR); ITERS := StringToInteger(ITERS);

A := QuaternionAlgebra(D); O := QuaternionOrder(A, N);
G := FuchsianGroup(O);
H := UpperHalfPlane();
CC := ComplexField(PREC);

ctr := H ! (9/10*H.1);          // Magma's default
printf "=== D=%o N=%o vol %o ===\n", D, N, RealField(5)!ArithmeticVolume(G);
best := -1; bestctr := ctr;
for it in [0..ITERS] do
    DDi := UnitDisc( : Center := ctr, Precision := PREC);
    t := Cputime();
    FDi := FundamentalDomain(G, DDi);
    tf := Cputime(t);
    vs := [ComplexValue(v) : v in FDi];
    mods := Sort([Abs(v) : v in vs]);
    rho := mods[#mods];
    printf "  step %o: centre %o, rho = %o (min vertex %o), %o sides, domain %o s\n",
        it, ComplexValue(ctr), RealField(6)!rho, RealField(4)!mods[1], #FDi, RealField(4)!tf;
    if best lt 0 or rho lt best then best := rho; bestctr := ctr; end if;
    if it eq ITERS then break; end if;
    c := (&+vs) / #vs;                       // Euclidean centroid of the vertices in the disc
    if Abs(c) lt 1e-12 then printf "  centroid already at the centre, stopping\n"; break; end if;
    ctr := DiscToPlane(H, DDi ! c);
end for;
printf "RECENTRE D=%o N=%o: best rho %o at centre %o (default gave the step-0 value)\n",
    D, N, RealField(6)!best, ComplexValue(bestctr);
quit;
