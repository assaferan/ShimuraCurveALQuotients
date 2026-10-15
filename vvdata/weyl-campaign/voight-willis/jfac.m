// Remark 3.3 of the paper: work with f_norm(z) = f(z) (Im z)^{k/2}, whose automorphy factor has
// ABSOLUTE VALUE 1, "since then the contribution from the automorphy factor may otherwise be quite
// large" when rho is large.  The current AutSystem uses the raw factor j^k.  Measure it.
AttachSpec("/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/ShimuraQuotients.spec");
SetColumns(0);
load "/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/vvdata/weyl-campaign/voight-willis/vw3.m";
load "/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/vvdata/weyl-campaign/voight-willis/vwal.m";
D := StringToInteger(DISC); N := StringToInteger(LEVEL); PREC := 80; k := 4;
CC := ComplexField(PREC); RR := RealField(PREC);
A := QuaternionAlgebra(D); O := QuaternionOrder(A, N); G := FuchsianGroup(O);
H := UpperHalfPlane();
_, DD, FD := RecentredDisc(G, H, PREC, 3);
rho := Max([Abs(ComplexValue(v)) : v in FD]);
p := ComplexValue(DiscToPlane(H, DD ! 0) : Precision := PREC);
mrep := FuchsianMatrixRepresentation(A : Precision := PREC);
samp := ReducedSamples(G, O, DD, mrep, p, rho, 200, CC);
js := Sort([Abs(DiscJ(s[3], s[1]))^k : s in samp]);
med := func<L | L[(#L+1) div 2]>;
printf "D=%o N=%o rho %o k=%o samples %o: |j|^k  min %o  med %o  MAX %o   (spread %o)\n",
   D, N, RR!rho, k, #samp, RR!js[1], RR!med(js), RR!js[#js], RR!(js[#js]/js[1]);
printf "JFAC DONE\n";
quit;
