// Calibrate the Hecke relation of eq. (18) against a known eigenvalue before trusting it anywhere.
//
// At D = 6, N = 1, weight 4 the space is one-dimensional, so the expansion the Gamma relations
// already produce IS an eigenform, and the ratio
//
//      ( sum_n K_n^h b_n ) / b_0
//
// is the eigenvalue in whatever normalisation eq. (18) carries.  Jacquet-Langlands identifies the
// space with the 6-new part of S_4(Gamma_0(6)), which is one-dimensional, so Magma's modular
// symbols give the classical a_p to compare against.  The ratio of the two is the normalisation,
// and it must be the SAME power of p at every prime or the construction is wrong.
//
//   magma -b DISC:=6 LEVEL:=1 PR:=150 TGT:=20 QMUL:=2 K:=4 RC:=3 PS:=5,7,11,13 heckecal.m < /dev/null
AttachSpec("/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/ShimuraQuotients.spec");
SetColumns(0);
load "/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/vvdata/weyl-campaign/voight-willis/vw3.m";
load "/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/vvdata/weyl-campaign/voight-willis/vwal.m";
load "/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/vvdata/weyl-campaign/voight-willis/vwhecke.m";

D := StringToInteger(DISC); N := StringToInteger(LEVEL);
PREC := StringToInteger(PR); TARGET := StringToInteger(TGT);
k := StringToInteger(K); rcit := StringToInteger(RC);
qmul := StringToInteger(QMUL);
ps := [StringToInteger(x) : x in Split(PS, ",")];
CC := ComplexField(PREC);
R8 := RealField(8);
R4 := RealField(4);

A := QuaternionAlgebra(D); O := QuaternionOrder(A, N);
G := FuchsianGroup(O);
H := UpperHalfPlane();
if rcit gt 0 then _, DD, FD := RecentredDisc(G, H, PREC, rcit);
else DD := UnitDisc( : Precision := PREC); FD := FundamentalDomain(G, DD); end if;
rho := Max([Abs(ComplexValue(v)) : v in FD]);
p0 := ComplexValue(DiscToPlane(H, DD ! 0) : Precision := PREC);
mrep := FuchsianMatrixRepresentation(A : Precision := PREC);

X1 := CreateShimuraQuot(D, N, {Integers() | 1});
g1, es1 := QuotientSignature(X1);
expd := DimMkSig(k, g1, es1);

NT := TermsFor(rho, TARGET);
Q := qmul*NT;
samp := ReducedSamples(G, O, DD, mrep, p0, rho, Q, CC);
if #samp lt qmul*NT and #samp gt 0 then
    Q := Ceiling(1.15*qmul*NT*Q/#samp);
    samp := ReducedSamples(G, O, DD, mrep, p0, rho, Q, CC);
end if;
Raw := AutSystem(samp, k, NT, CC);
nms := [ Sqrt(&+[ Abs(Raw[i][j])^2 : i in [1..#samp] ]) : j in [1..NT] ];
bs, smin := SmallestVec(NormaliseColumns(Raw));

printf "=== D=%o N=%o k=%o rho %o NT %o rows %o prec %o ===\n", D, N, k, R4!rho, NT, #samp, PREC;
printf "signature (%o; %o), expected dim %o, smallest singular value %o\n", g1, es1, expd, R4!smin;
if expd ne 1 then
    printf "WARNING: the space is %o-dimensional, so the kernel vector is not an eigenform and\n", expd;
    printf "         the ratio below is not an eigenvalue.  Calibrate at a one-dimensional base.\n";
end if;

// the raw system's kernel coefficients
braw := [ CC | bs[j]/nms[j] : j in [1..NT] ];
sc := braw[1];
braw := [ b/sc : b in braw ];           // normalise b_0 = 1
printf "b_0 normalised to 1; |b_1| %o |b_2| %o |b_%o| %o\n",
    R4!Abs(braw[2]), R4!Abs(braw[3]), NT-1, R4!Abs(braw[NT]);

// residual of the Gamma relations on this vector, as a sanity check that it really is the form
res := Max([ Abs(&+[ (samp[i][2]^n - DiscJ(samp[i][3], samp[i][1])^k * samp[i][1]^n) * braw[n+1]
                     : n in [0..NT-1] ]) : i in [1..#samp] ]);
printf "max residual of the automorphy relations on b: %o\n", R4!res;

// Jacquet-Langlands: M_k(Gamma) here is the D-new part of S_k(Gamma_0(DN)).  Calibration is only
// meaningful where that is one-dimensional, so insist on a single newform.
nfs := Newforms(CuspForms(Gamma0(N*D), k));
printf "S_%o(Gamma_0(%o)) has %o newform orbit(s)\n", k, N*D, #nfs;
assert #nfs eq 1 and #nfs[1] eq 1;
fcl := nfs[1][1];
qex := qExpansion(fcl, Max(ps)+2);

for pp in ps do
    if D*N mod pp eq 0 then printf "  p=%o divides D*N, skipped\n", pp; continue; end if;
    reps := HeckeReps(O, pp);
    if #reps ne pp+1 then
        printf "  p=%o: found %o coset representatives, expected %o -- widen the box\n",
            pp, #reps, pp+1;
        continue;
    end if;
    Kv, tgts := HeckeKernelRow(G, O, DD, mrep, p0, reps, k, NT, CC);
    rat := HeckeRatio(Kv, braw);
    apc := Integers() ! Coefficient(qex, pp);
    maxw := Max([ Abs(t[1]) : t in tgts ]);
    maxr := Max([ Abs(t[2]) : t in tgts ]);
    printf "  p=%o: %o cosets, |Hecke image| max %o -> reduced max %o (rho %o)\n",
        pp, #reps, R4!maxw, R4!maxr, R4!rho;
    printf "        ratio %o   classical a_%o = %o   ratio/a_p = %o   ratio/(a_p*p^(k/2-1)) = %o\n",
        R8!Real(rat), pp, apc, R8!(Real(rat)/apc), R8!(Real(rat)/(apc*pp^(k div 2 - 1)));
    printf "        |Im(ratio)| %o\n", R4!Abs(Imaginary(rat));
end for;
printf "HECKECAL DONE\n";
quit;
