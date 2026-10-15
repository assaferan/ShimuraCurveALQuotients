// Validation of the Atkin-Lehner layer: for every Atkin-Lehner subgroup W of a base, the dimension
// of the W-invariant space of weight-k forms computed from the expansions, against the dimension
// predicted from the quotient's own signature (genus from the genus formula, elliptic points from
// NumberOfEllipticPointsByCMOrder).  At k = 2 the prediction is the genus of X_0^D(N)/W itself.
//
//   magma -b DISC:=6 LEVEL:=1 PR:=100 NTERMS:=50 QQ:=60 BOX:=3 aldims.m < /dev/null
AttachSpec("ShimuraQuotients.spec");
SetColumns(0);
load "/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/vvdata/weyl-campaign/voight-willis/vw3.m";
load "/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/vvdata/weyl-campaign/voight-willis/vwal.m";

D := StringToInteger(DISC); N := StringToInteger(LEVEL);
PREC := StringToInteger(PR);
BOXSZ := StringToInteger(BOX); TARGET := StringToInteger(TGT);
DN := D*N;
CC := ComplexField(PREC);

t0 := Cputime();
A := QuaternionAlgebra(D); O := QuaternionOrder(A, N);
G := FuchsianGroup(O);
H := UpperHalfPlane();
_, DD, FD := RecentredDisc(G, H, PREC, 3);
rho := Max([Abs(ComplexValue(v)) : v in FD]);
p := ComplexValue(DiscToPlane(H, DD ! 0) : Precision := PREC);
mrep := FuchsianMatrixRepresentation(A : Precision := PREC);
tdom := Cputime(t0);

// terms and tolerance follow the domain radius, not a fixed choice
NT := TermsFor(rho, TARGET);
tol := ZeroLevel(rho, NT : prec := PREC);
Q := 3*NT;
sampG := ReducedSamples(G, O, DD, mrep, p, rho, Q, CC);
printf "D=%o N=%o DN=%o vol %o rho=%o domain %o s, %o sides, N=%o tol=%o, %o Gamma samples, prec %o\n",
    D, N, DN, RealField(5)!ArithmeticVolume(G), RealField(6)!rho, RealField(4)!tdom, #FD, NT,
    RealField(3)!tol, #sampG, PREC;

// the Atkin-Lehner elements, once per base
mu := AssociativeArray();
for m in [d : d in Divisors(DN) | d gt 1 and GCD(d, DN div d) eq 1] do
    ok, x := ALElement(O, m, BOXSZ);
    error if not ok, Sprintf("no Atkin-Lehner element of reduced norm %o found in a box of %o", m, BOXSZ);
    mu[m] := x;
end for;
printf "Atkin-Lehner elements found for m = %o\n", Sort(SetToSequence(Keys(mu)));

nok := 0; ntot := 0;
for Wp in ALSubgroups(DN) do
    W := Wp[1];
    X := CreateShimuraQuot(D, N, W);
    g, es := QuotientSignature(X);
    gens := ALGens(W, DN);
    samp := sampG cat (#gens eq 0 select [] else ALSamples(G, O, DD, mrep, p, rho, Q, CC, [mu[m] : m in gens]));
    line := Sprintf("  W=%-26o g=%o sig(%o;%o) samples %o |", Sort(SetToSequence(W)), g, g, es, #samp);
    for k in [2, 4, 6] do
        exp := DimMkSig(k, g, es);
        sv := SingVals(NormaliseColumns(AutSystem(samp, k, NT, CC)));
        got, gap := GapDim(sv, tol);
        ntot +:= 1; if got eq exp then nok +:= 1; end if;
        line cat:= Sprintf(" k=%o exp %o got %o%o (gap %o)", k, exp, got,
            got eq exp select "" else " MISMATCH", RealField(3)!gap);
    end for;
    print line;
end for;
printf "ALDIMS D=%o N=%o: %o of %o dimensions correct\n", D, N, nok, ntot;
quit;
