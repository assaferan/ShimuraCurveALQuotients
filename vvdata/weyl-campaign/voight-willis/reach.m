// Does a base resolve, and what does it cost?  Two knobs, measured separately, following the
// diagnosis of 2026-10-15:
//
//   RADF   the sampling circle, as a percentage of the domain's circumradius rho.  Points outside
//          the domain are the ones that carry a relation, and at RADF = 100 essentially the whole
//          circle is outside once the reduction works, so going above 100 buys nothing and costs
//          terms: the series is evaluated out to the SAMPLING radius, so the term count follows
//          that and not rho, and it grows without bound as the radius approaches 1.
//
//   TGTS   a list of truncation targets, so the term count can be sized by watching the gap grow
//          rather than by one fixed target.
//
// The diagnostic that tells the two failure modes apart is the smallest singular value against the
// truncation floor rad^NT, together with the size of the gap.  A spectrum sitting far BELOW the
// floor with no gap is a rank-deficient matrix, which more terms never repair; that is what the
// clustered sampling produced before the reduction was fixed.  A kernel a hundred times ABOVE the
// floor with a gap of 1e9 or more is a resolved base, the kernel being limited by the quality of
// the reduction rather than by truncation.
//
//   magma -b DISC:=15 LEVEL:=1 PR:=150 RADF:=110 QMUL:=3 K:=4 RC:=3 TGTS:=12,20,28 reach.m < /dev/null
AttachSpec("/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/ShimuraQuotients.spec");
SetColumns(0);
load "/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/vvdata/weyl-campaign/voight-willis/vw3.m";
load "/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/vvdata/weyl-campaign/voight-willis/vwal.m";

D := StringToInteger(DISC); N := StringToInteger(LEVEL);
PREC := StringToInteger(PR);
radf := StringToInteger(RADF)/100;
k := StringToInteger(K);
rcit := StringToInteger(RC);
qmul := StringToInteger(QMUL);
tgts := [StringToInteger(x) : x in Split(TGTS, ",")];
CC := ComplexField(PREC);
RR := RealField(PREC);
R4 := RealField(4);

A := QuaternionAlgebra(D); O := QuaternionOrder(A, N);
G := FuchsianGroup(O);
H := UpperHalfPlane();
t0 := Cputime();
if rcit gt 0 then _, DD, FD := RecentredDisc(G, H, PREC, rcit);
else DD := UnitDisc( : Precision := PREC); FD := FundamentalDomain(G, DD); end if;
tdom := Cputime(t0);
vmods := Sort([Abs(ComplexValue(v)) : v in FD]);
rho := vmods[#vmods];
rad := radf * rho;
p := ComplexValue(DiscToPlane(H, DD ! 0) : Precision := PREC);
mrep := FuchsianMatrixRepresentation(A : Precision := PREC);

X1 := CreateShimuraQuot(D, N, {Integers() | 1});
g1, es1 := QuotientSignature(X1);
expd := DimMkSig(k, g1, es1);

printf "=== D=%o N=%o vol %o sides %o (domain %o s) ===\n",
    D, N, R4!ArithmeticVolume(G), #FD, R4!tdom;
printf "signature (%o; %o), weight %o, expected dim %o\n", g1, es1, k, expd;
printf "rho %o, sampling radius %o (%o x rho), precision %o\n", R4!rho, R4!rad, R4!radf, PREC;

prev := 0.0;
for tgt in tgts do
    NT := TermsFor(rad, tgt);
    if NT ge PREC*2 then
        printf "  target 1e-%o: NT %o exceeds the working precision, skipped\n", tgt, NT;
        continue;
    end if;
    // Points that land inside the domain reduce to themselves and carry no relation, so the
    // usable rows are a fraction of Q that depends on the base.  Oversample until there are
    // qmul*NT rows rather than qmul*NT samples, measuring the retention on the first pass.
    Q := qmul*NT;
    t1 := Cputime();
    samp := ReducedSamples(G, O, DD, mrep, p, rad, Q, CC);
    ret := #samp/Q;
    if #samp lt qmul*NT and ret gt 0 then
        Q := Min(Ceiling(1.15*qmul*NT/ret), 40*NT);
        samp := ReducedSamples(G, O, DD, mrep, p, rad, Q, CC);
    end if;
    tsamp := Cputime(t1);
    if #samp lt NT then
        printf "  target 1e-%o: NT %o but only %o rows from %o samples (retention %o%%), underdetermined\n",
            tgt, NT, #samp, Q, R4!(100*ret);
        continue;
    end if;
    t2 := Cputime();
    sv := SingVals(NormaliseColumns(AutSystem(samp, k, NT, CC)));
    tsvd := Cputime(t2);
    floor := rad^NT;
    gr, grat := GapDimStrict(sv, floor, 1000);
    cond := sv[#sv]/sv[1];
    printf "  1e-%o: NT %o rows %o (%o%% of %o) | floor %o smallest %o (x%o) | cond %o | dim %o%o gap %o | %o+%o s\n",
        tgt, NT, #samp, Round(100*#samp/Q), Q, R4!floor, R4!sv[1], R4!(sv[1]/floor), R4!cond,
        gr lt 0 select "UNRESOLVED" else Sprint(gr), gr eq expd select "" else " MISMATCH", R4!grat, R4!tsamp, R4!tsvd;
    prev := grat;
end for;
printf "REACH DONE\n";
quit;
