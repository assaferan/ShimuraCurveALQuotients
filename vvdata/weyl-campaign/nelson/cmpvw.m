// Stage 1: compare Nelson's evaluator against the Voight-Willis expansions at the same points.
//
// At X_0^D(1) with the weight chosen so that the space of holomorphic forms is ONE dimensional,
// any two nonzero forms are proportional, so the test is whether the RATIO f(z_2)/f(z_1) agrees.
// Three points would already over-determine nothing; every point past the first is a prediction.
//
//   magma -b DISC:=6 K:=4 PR:=30 CUT:=35 TGT:=20 NPTS:=8 cmpvw.m < /dev/null
AttachSpec("/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/ShimuraQuotients.spec");
SetColumns(0);
VW := "/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/vvdata/weyl-campaign/voight-willis/";
load "/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/vvdata/weyl-campaign/voight-willis/vw3.m";
load "/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/vvdata/weyl-campaign/voight-willis/vwal.m";
load "/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/vvdata/weyl-campaign/nelson/nelson.m";

D := StringToInteger(DISC); N := 1; k := StringToInteger(K);
PREC := StringToInteger(PR); cutv := StringToInteger(CUT);
TARGET := StringToInteger(TGT); npts := StringToInteger(NPTS);
CC := ComplexField(PREC); RR := RealField(PREC); Ic := CC.1;

A := QuaternionAlgebra(D); O := QuaternionOrder(A, N);
G := FuchsianGroup(O);
H := UpperHalfPlane();
_, DD, FD := RecentredDisc(G, H, PREC, 3);
p := ComplexValue(DiscToPlane(H, DD ! 0) : Precision := PREC);
mrep := FuchsianMatrixRepresentation(A : Precision := PREC);
rho := Max([Abs(ComplexValue(v)) : v in FD]);
NT := TermsFor(rho, TARGET);
printf "=== D=%o k=%o  rho=%o  NT=%o  prec=%o  cut=%o\n", D, k, RR!rho, NT, PREC, cutv;

// ---- Voight-Willis side: the one-dimensional space of weight-k forms
t0 := Cputime();
samp := ReducedSamples(G, O, DD, mrep, p, rho, 3*NT, CC);
Mx := AutSystem(samp, k, NT, CC);
cn := [ Sqrt(&+[ Abs(Mx[i][j])^2 : i in [1..Nrows(Mx)] ]) : j in [1..NT] ];
v, smallest := SmallestVec(NormaliseColumns(Mx));
bs := [ v[j]/cn[j] : j in [1..NT] ];
sv := SingVals(NormaliseColumns(Mx));
printf "VW: %o samples, smallest sv %o, next %o, gap %o  (%o s)\n",
       #samp, RR!sv[1], RR!sv[2], RR!(sv[2]/sv[1]), RR!Cputime(t0);

// ---- Nelson side
f := Newforms(CuspForms(Gamma0(D*N), k))[1][1];
nmax := 6000;
qq := qExpansion(f, nmax+1);
afs := [ Rationals() ! Coefficient(qq, n) : n in [1..nmax] ];
af := func< n | n le nmax select CC ! afs[n] else 0 >;
printf "newform a_f(1..8): %o\n", afs[1..8];
dat := NelsonData(D, N, k, PREC : OO := O, MREP := mrep);

// ---- the sample points: disc points of modest modulus, pulled back to H
ws := [ CC | 0.42*rho*Exp(2*Pi(RR)*Ic*j/npts) + 0.08*rho : j in [0..npts-1] ];
zs := [ ComplexValue(DiscToPlane(H, DD ! w) : Precision := PREC) : w in ws ];

S11 := NelsonPair(dat, zs[1], zs[1], af, RR!cutv, CC);
F1 := EvalForm(bs, k, ws[1]);
printf "\n  j   |w_j|      |Nelson ratio|            VW ratio agreement\n";
errs := [];
for j in [2..npts] do
    Sj := NelsonPair(dat, zs[1], zs[j], af, RR!cutv, CC);
    rN := (Imaginary(zs[1])/Imaginary(zs[j]))^(k/2) * Sj / S11;     // f(z_j)/f(z_1), Nelson
    rV := EvalForm(bs, k, ws[j]) / F1;                               // f(z_j)/f(z_1), Voight-Willis
    e := Abs(rN - rV)/Abs(rV);
    Append(~errs, e);
    printf "  %-3o %-10o %-24o rel.diff %o\n", j, RR!Abs(ws[j]), CC!rN, RR!e;
end for;
printf "\nworst relative difference over %o predictions: %o\n", npts-1, RR!Max(errs);
printf "CMPVW DONE\n";
quit;
