// Why does resolution fail as the domain radius grows?  Separates the two candidate causes:
//
//   TRUNCATION -- the series is cut at NT terms, so nothing below the floor rho^NT is meaningful
//     and a kernel vector cannot be told from noise.  Signature: as NT and the working precision
//     rise together, the whole spectrum slides DOWN with the floor and a gap eventually appears.
//
//   CONDITIONING -- the rows are nearly dependent regardless of the floor, so there is a continuum
//     of near-null directions.  Signature: the spectrum does NOT follow the floor; the smallest
//     singular values stall far above rho^NT and no gap ever opens.
//
// Also reports how far the reduction actually moves a sample point, since a relation
// f(delta w) = j^k f(w) with delta w close to w carries almost no information however many rows
// there are.
//
//   magma -b DISC:=22 LEVEL:=1 PR:=120 TGT:=12 QMUL:=3 K:=4 RC:=3 walldiag.m < /dev/null
AttachSpec("/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/ShimuraQuotients.spec");
SetColumns(0);
load "/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/vvdata/weyl-campaign/voight-willis/vw3.m";
load "/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/vvdata/weyl-campaign/voight-willis/vwal.m";

D := StringToInteger(DISC); N := StringToInteger(LEVEL);
PREC := StringToInteger(PR); TARGET := StringToInteger(TGT);
k := StringToInteger(K); rcit := StringToInteger(RC);
CC := ComplexField(PREC);
RR := RealField(PREC);

A := QuaternionAlgebra(D); O := QuaternionOrder(A, N);
G := FuchsianGroup(O);
H := UpperHalfPlane();
if rcit gt 0 then _, DD, FD := RecentredDisc(G, H, PREC, rcit);
else DD := UnitDisc( : Precision := PREC); FD := FundamentalDomain(G, DD); end if;
vmods := Sort([Abs(ComplexValue(v)) : v in FD]);
rho := vmods[#vmods];
rad := rho;
p := ComplexValue(DiscToPlane(H, DD ! 0) : Precision := PREC);
mrep := FuchsianMatrixRepresentation(A : Precision := PREC);

NT := TermsFor(rad, TARGET);
Q := StringToInteger(QMUL)*NT;
floor := rho^NT;

t0 := Cputime();
samp := ReducedSamples(G, O, DD, mrep, p, rad, Q, CC);
tsamp := Cputime(t0);

// how far does the reduction move the point?  |w| is rad for every sample.
movs := Sort([Abs(s[2] - s[1]) : s in samp]);
tgts := Sort([Abs(s[2]) : s in samp]);
med := func<L | L[(#L+1) div 2]>;

X1 := CreateShimuraQuot(D, N, {Integers() | 1});
g1, es1 := QuotientSignature(X1);
expd := DimMkSig(k, g1, es1);

t1 := Cputime();
Mx := NormaliseColumns(AutSystem(samp, k, NT, CC));
sv := SingVals(Mx);
tsvd := Cputime(t1);

cond := sv[#sv]/sv[1];
gd, gap := GapDim(sv, ZeroLevel(rad, NT : prec := PREC));
gr, grat := GapDimRel(sv, 1e-3);

printf "D=%o N=%o k=%o | rho %o rad %o NT %o prec %o rows %o\n",
   D, N, k, RR!rho, RR!rad, NT, PREC, #samp;
printf "  truncation floor rho^NT = %o\n", RR!floor;
printf "  smallest sv = %o   (sv_min / floor = %o)\n", RR!sv[1], RR!(sv[1]/floor);
printf "  largest sv  = %o   condition number = %o\n", RR!sv[#sv], RR!cond;
printf "  10 smallest sv: %o\n", [RR!s : s in sv[1..Min(10,#sv)]];
printf "  expected dim %o | GapDim %o (gap %o) | GapDimRel %o (gap %o)\n",
   expd, gd, RR!gap, gr, RR!grat;
printf "  reduction displacement |w'-w|: min %o med %o max %o\n", RR!movs[1], RR!med(movs), RR!movs[#movs];
printf "  target moduli |w'|:            min %o med %o max %o   (rad = %o)\n",
   RR!tgts[1], RR!med(tgts), RR!tgts[#tgts], RR!rad;
printf "  times: samples %o s, svd %o s\n", RR!tsamp, RR!tsvd;
printf "WALLDIAG DONE\n";
quit;
