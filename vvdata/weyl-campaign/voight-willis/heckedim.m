// Does the Hecke relation of eq. (18) rescue a base where the automorphy relations alone do not
// resolve?  This is the paper's main lever.  A d-dimensional space asks the linear system for a
// d-dimensional kernel separated from the rest of the spectrum; one Hecke relation per prime,
// with the eigenvalue of a chosen eigenform substituted, asks instead for d separate
// ONE-dimensional kernels, and a one-dimensional kernel survives far worse conditioning.
//
// The eigenvalues come from the classical side: Jacquet-Langlands identifies M_k(Gamma) with the
// D-new part of S_k(Gamma_0(DN)), so Magma's modular symbols supply them, and each complex
// embedding of the coefficient field gives one eigenform on the quaternion side.  The
// normalisation is fixed by heckecal.m: the operator of (17), with the coset matrices scaled to
// determinant 1, has eigenvalue p^(1-k/2) a_p.
//
//   magma -b DISC:=15 LEVEL:=1 PR:=200 TGT:=20 RADF:=100 QMUL:=2 K:=4 RC:=3 PS:=2,7,11 heckedim.m < /dev/null
AttachSpec("/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/ShimuraQuotients.spec");
SetColumns(0);
load "/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/vvdata/weyl-campaign/voight-willis/vw3.m";
load "/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/vvdata/weyl-campaign/voight-willis/vwal.m";
load "/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/vvdata/weyl-campaign/voight-willis/vwhecke.m";

D := StringToInteger(DISC); N := StringToInteger(LEVEL);
PREC := StringToInteger(PR); TARGET := StringToInteger(TGT);
radf := StringToInteger(RADF)/100;
k := StringToInteger(K); rcit := StringToInteger(RC);
qmul := StringToInteger(QMUL);
ps := [StringToInteger(x) : x in Split(PS, ",")];
CC := ComplexField(PREC);
R4 := RealField(4);

A := QuaternionAlgebra(D); O := QuaternionOrder(A, N);
G := FuchsianGroup(O); H := UpperHalfPlane();
if rcit gt 0 then _, DD, FD := RecentredDisc(G, H, PREC, rcit);
else DD := UnitDisc( : Precision := PREC); FD := FundamentalDomain(G, DD); end if;
rho := Max([Abs(ComplexValue(v)) : v in FD]);
rad := radf*rho;
p0 := ComplexValue(DiscToPlane(H, DD ! 0) : Precision := PREC);
mrep := FuchsianMatrixRepresentation(A : Precision := PREC);

X1 := CreateShimuraQuot(D, N, {Integers() | 1});
g1, es1 := QuotientSignature(X1);
expd := DimMkSig(k, g1, es1);

NT := TermsFor(rad, TARGET);
Q := qmul*NT;
t0 := Cputime();
samp := ReducedSamples(G, O, DD, mrep, p0, rad, Q, CC);
if #samp lt qmul*NT and #samp gt 0 then
    Q := Min(Ceiling(1.15*qmul*NT*Q/#samp), 40*NT);
    samp := ReducedSamples(G, O, DD, mrep, p0, rad, Q, CC);
end if;
tsamp := Cputime(t0);

printf "=== D=%o N=%o k=%o | rho %o rad %o (%ox) NT %o rows %o of %o prec %o ===\n",
    D, N, k, R4!rho, R4!rad, R4!radf, NT, #samp, Q, PREC;
printf "signature (%o; %o), expected dim %o  (sampling %o s)\n", g1, es1, expd, R4!tsamp;

Raw := AutSystem(samp, k, NT, CC);
floor := rad^NT;
t1 := Cputime();
sv0 := SingVals(NormaliseColumns(Raw));
d0, gap0 := GapDimStrict(sv0, floor, 1000);
printf "Gamma relations alone: floor %o, smallest sv %o (x%o), dim %o%o, gap %o  (%o s)\n",
    R4!floor, R4!sv0[1], R4!(sv0[1]/floor),
    d0 lt 0 select "UNRESOLVED (no gap)" else Sprint(d0),
    d0 eq expd select "" else " MISMATCH", R4!gap0, R4!Cputime(t1);
printf "  6 smallest: %o\n", [R4!s : s in sv0[1..Min(6,#sv0)]];

// classical eigenvalues, one eigenform per complex embedding of each newform's coefficient field
Snew := NewSubspace(CuspidalSubspace(ModularSymbols(N*D, k, 1)), D);
dec := NewformDecomposition(Snew);
forms := [];
for S in dec do
    f := qEigenform(S, Max(ps)+1);
    Kf := BaseRing(Parent(f));
    if Kf cmpeq Rationals() then
        Append(~forms, [ CC | Rationals() ! Coefficient(f, pp) : pp in ps ]);
    else
        for r in [ z[1] : z in Roots(DefiningPolynomial(Kf), CC) ] do
            Append(~forms, [ CC | &+[ (CC ! Eltseq(Coefficient(f, pp))[i]) * r^(i-1)
                                      : i in [1..Degree(Kf)] ] : pp in ps ]);
        end for;
    end if;
end for;
printf "classical D-new space: %o newform orbit(s), %o eigenforms, expected %o\n",
    #dec, #forms, expd;

// one Hecke row per prime, built once and reused for every eigenform
Ks := [];
for pp in ps do
    if (D*N) mod pp eq 0 then printf "  (p=%o divides D*N, skipped)\n", pp; continue; end if;
    reps := HeckeReps(O, pp);
    if #reps ne pp+1 then
        printf "  (p=%o: %o cosets found, expected %o, skipped)\n", pp, #reps, pp+1; continue;
    end if;
    Kv, tg := HeckeKernelRow(G, O, DD, mrep, p0, reps, k, NT, CC);
    Append(~Ks, <pp, Kv, Max([Abs(t[2]) : t in tg])>);
    printf "  p=%o: %o cosets, worst reduced Hecke image %o (rho %o), its truncation %o\n",
        pp, #reps, R4!Ks[#Ks][3], R4!rho, R4!(Ks[#Ks][3]^NT);
end for;

nms := [ Sqrt(&+[ Abs(Raw[i][j])^2 : i in [1..#samp] ]) : j in [1..NT] ];
Gm := NormaliseColumns(Raw);

// shift = 0 is the true eigenvalue; a nonzero shift is a negative control, and the dimension
// must drop to 0 there or the Hecke rows are not constraining anything.
function WithHecke(Gm, Ks, ps, evs, k, NT, nms, shift, CC)
    rows := [];
    for K in Ks do
        pp := K[1];
        ev := evs[Index(ps, pp)] * (CC ! pp)^(1 - k div 2) + shift;
        rel := HeckeRelation(K[2], ev);
        r := [ rel[j]/nms[j] : j in [1..NT] ];
        nr := Sqrt(&+[ Abs(x)^2 : x in r ]);
        if nr gt 0 then r := [ x/nr : x in r ]; end if;   // one equation, rescaled to match the
        Append(~rows, r);                                  // Gamma rows in weight
    end for;
    return SingVals(VerticalJoin(Gm, Matrix(CC, #rows, NT, &cat rows)));
end function;

for fi in [1..#forms] do
    if #Ks eq 0 then break; end if;
    sv := WithHecke(Gm, Ks, ps, forms[fi], k, NT, nms, CC ! 0, CC);
    dd, gp := GapDimStrict(sv, floor, 1000);
    svc := WithHecke(Gm, Ks, ps, forms[fi], k, NT, nms, CC ! 1, CC);
    dc := GapDimStrict(svc, floor, 1000);
    printf "eigenform %o (a_%o = %o): dim %o%o, smallest sv %o, gap %o | control (a_p+1) dim %o%o\n",
        fi, ps[1], R4!Real(forms[fi][1]),
        dd lt 0 select "UNRESOLVED (no gap)" else Sprint(dd),
        dd eq 1 select "" else " NOT ONE-DIMENSIONAL",
        R4!sv[1], R4!gp, dc lt 0 select "none" else Sprint(dc),
        dc le 0 select "" else " CONTROL FAILED";
    printf "   6 smallest: %o\n", [R4!s : s in sv[1..Min(6,#sv)]];
    printf "   control    : %o   (smallest sv rises by x%o)\n",
        [R4!s : s in svc[1..Min(6,#svc)]], R4!(svc[1]/sv[1]);
end for;
printf "HECKEDIM DONE\n";
quit;
