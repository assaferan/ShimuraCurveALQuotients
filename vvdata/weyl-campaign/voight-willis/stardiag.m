AttachSpec("/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/ShimuraQuotients.spec");
SetColumns(0);
load "/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/vvdata/weyl-campaign/voight-willis/vw3.m";
load "/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/vvdata/weyl-campaign/voight-willis/vwal.m";
load "/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/vvdata/weyl-campaign/voight-willis/vwmodel.m";

D := StringToInteger(DISC); N := StringToInteger(LEVEL);
PREC := StringToInteger(PR); TGT := StringToInteger(T);
CC := ComplexField(PREC); RR := RealField(PREC); R4 := RealField(4);
DN := D*N;
A := QuaternionAlgebra(D); O := QuaternionOrder(A, N);
G := FuchsianGroup(O); H := UpperHalfPlane();
_, DD, FD := RecentredDisc(G, H, PREC, 3);
rhoG := Max([Abs(ComplexValue(v)) : v in FD]);
pt := ComplexValue(DiscToPlane(H, DD ! 0) : Precision := PREC);
mrep := FuchsianMatrixRepresentation(A : Precision := PREC);
Ws := SidePairingMatrices(G, DD, mrep, pt, CC);
Wstar := { Integers() | d : d in Divisors(DN) | GCD(d, DN div d) eq 1 };
Als := [];
for m in Sort([ d : d in Wstar | d gt 1 ]) do
    ok, x := ALElement(O, m, 3);
    Append(~Als, <m, DiscMatrix(mrep, x, CC ! pt, CC)>);
end for;
Xs := CreateShimuraQuot(D, N, Wstar); gs, ess := QuotientSignature(Xs);
rmax := RR ! 0;
for i in [1..2000] do
    r := (RR ! (0.05 + 0.95*((i*5) mod 23)/23)) * rhoG;
    w := CC ! (r * Exp(2*Pi(RR)*CC.1*i*(Sqrt(RR ! 5)-1)/2));
    _, v := ReduceOverStar(Ws, Als, w, CC);
    if Abs(v) gt rmax then rmax := Abs(v); end if;
end for;
printf "D=%o rho(Gamma)=%o rho(Gamma*)=%o star sig (%o;%o)\n", D, R4!rhoG, R4!rmax, gs, ess;

rad := Min(RR ! 0.999, 1.03*rmax);
NT := TermsFor(rad, TGT);
samp := StarSamples(Ws, Als, rad, 30*NT, CC);
for am in Als do samp cat:= StarSamplesVia(Ws, Als, am[2], am[1], rad, 6*NT, CC); end for;
// how many DISTINCT group elements do the reductions realise?
els := [];
for s in samp do
    if exists{ M : M in els | Abs(M[1][1]-s[3][1][1]) + Abs(M[1][2]-s[3][1][2]) lt 1e-20 } then continue; end if;
    Append(~els, s[3]);
    if #els gt 400 then break; end if;
end for;
printf "NT %o, %o samples, at least %o distinct group elements realised\n", NT, #samp, #els;
sgnS := AssociativeArray();
for m in Wstar do sgnS[m] := 1; end for;
k := 4; d := DimMkSig(k, gs, ess);
for mult in [2, 6, 16] do
    M, nr, jlo, jhi := SystemLabelled(samp, k, NT, sgnS, mult*NT, CC);
    sv := SingVals(NormaliseColumns(M));
    printf "  rows %o (%ox NT) floor %o | k=%o expect %o : %o\n",
        nr, mult, R4!(rad^NT), k, d, [R4!s : s in sv[1..Min(d+3, #sv)]];
end for;
// control: the GAMMA-only system on the Gamma circle at the same target, old style
sampG := ReducedSamples(G, O, DD, mrep, pt, rhoG, 6*TermsFor(rhoG, TGT), CC);
NTg := TermsFor(rhoG, TGT);
X1 := CreateShimuraQuot(D, N, {Integers()|1}); g1, es1 := QuotientSignature(X1);
for mult in [2, 6] do
    rows := [];
    for s in sampG do
        w := s[1]; wp := s[2]; jk := DiscJ(s[3], w)^k;
        Append(~rows, [ wp^n - jk*w^n : n in [0..NTg-1] ]);
        if #rows ge mult*NTg then break; end if;
    end for;
    Mg := Matrix(CC, #rows, NTg, &cat rows);
    sv := SingVals(NormaliseColumns(Mg));
    printf "  CONTROL Gamma-only NTg %o rows %o floor %o | k=%o expect %o : %o\n",
        NTg, #rows, R4!(rhoG^NTg), k, DimMkSig(k, g1, es1),
        [R4!s : s in sv[1..Min(DimMkSig(k,g1,es1)+3, #sv)]];
end for;
printf "STARDBG2 DONE\n";
quit;
