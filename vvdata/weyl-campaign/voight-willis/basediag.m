// Does the BASE library (Gamma relations only, no Atkin-Lehner) resolve at a given base?
// Separates "the expansion library does not scale past D = 6" from "the Atkin-Lehner layer is
// wrong".  Also reports the sampling radius actually used and the moduli of the reduced targets.
//
//   magma -b DISC:=10 LEVEL:=1 PR:=120 TGT:=12 RADF:=100 RC:=3 QMUL:=3 KS:=2,4 basediag.m < /dev/null
//   (RADF in percent of the maximal vertex modulus; RC = recentring iterations, 0 for Magma's default)
AttachSpec("ShimuraQuotients.spec");
SetColumns(0);
load "/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/vvdata/weyl-campaign/voight-willis/vw3.m";
load "/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/vvdata/weyl-campaign/voight-willis/vwal.m";

D := StringToInteger(DISC); N := StringToInteger(LEVEL);
PREC := StringToInteger(PR); TARGET := StringToInteger(TGT);
radf := StringToInteger(RADF)/100;   // percentage of the maximal vertex modulus
ks := [StringToInteger(x) : x in Split(KS, ",")];
CC := ComplexField(PREC);
tiny := RealField(PREC) ! 10^-10;

A := QuaternionAlgebra(D); O := QuaternionOrder(A, N);
G := FuchsianGroup(O);
H0 := UpperHalfPlane();
rcit := StringToInteger(RC);
if rcit gt 0 then
    _, DD, FD := RecentredDisc(G, H0, PREC, rcit);
else
    DD := UnitDisc( : Precision := PREC);
    FD := FundamentalDomain(G, DD);
end if;
vmods := Sort([Abs(ComplexValue(v)) : v in FD]);
rhomax := vmods[#vmods];
rad := radf * rhomax;

H := UpperHalfPlane();
p := ComplexValue(DiscToPlane(H, DD ! 0) : Precision := PREC);
mrep := FuchsianMatrixRepresentation(A : Precision := PREC);

// terms sized to the SAMPLING radius actually used
NT := TermsFor(rad, TARGET);
tol := ZeroLevel(rad, NT : prec := PREC);
Q := StringToInteger(QMUL)*NT;

printf "=== D=%o N=%o vol %o sides %o ===\n", D, N, RealField(5)!ArithmeticVolume(G), #FD;
printf "vertex moduli: min %o max %o   sampling radius %o (= %o * rhomax)\n",
    RealField(4)!vmods[1], RealField(4)!rhomax, RealField(5)!rad, RealField(3)!radf;
printf "NT=%o tol=%o Q=%o prec=%o\n", NT, RealField(3)!tol, Q, PREC;

t0 := Cputime();
samp := ReducedSamples(G, O, DD, mrep, p, rad, Q, CC);
tms := Sort([Abs(s[2]) : s in samp]);
nover := #[x : x in tms | x gt rhomax + tiny];
printf "Gamma samples %o in %o s; target moduli min %o med %o max %o; over rhomax: %o\n",
    #samp, RealField(4)!Cputime(t0), RealField(4)!tms[1],
    RealField(4)!tms[(#tms+1) div 2], RealField(4)!tms[#tms], nover;

X1 := CreateShimuraQuot(D, N, {Integers() | 1});
g1, es1 := QuotientSignature(X1);
printf "signature (%o; %o)\n", g1, es1;
for k in ks do
    exp := DimMkSig(k, g1, es1);
    t1 := Cputime();
    sv := SingVals(NormaliseColumns(AutSystem(samp, k, NT, CC)));
    got, gap := GapDim(sv, tol);
    printf "  k=%o exp %o got %o%o gap %o (%o s) | 8 smallest %o\n", k, exp, got,
        got eq exp select "" else "  MISMATCH", RealField(3)!gap, RealField(4)!Cputime(t1),
        [RealField(3)!s : s in sv[1..Min(12, #sv)]];
end for;
printf "BASEDIAG DONE\n";
quit;
