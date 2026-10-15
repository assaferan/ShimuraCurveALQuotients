// Diagnostic for the D = 10 failure of the Atkin-Lehner layer.  Three questions, in order:
//   (1) does the BASE library (trivial W, Gamma relations only) resolve at this base at all?
//   (2) what is the distribution of |w| of the points actually used, for Gamma and for AL samples?
//   (3) how does the singular-value gap respond when rows with a large reduced modulus are dropped?
//
//   magma -b DISC:=10 LEVEL:=1 PR:=150 TGT:=14 BOX:=3 aldiag.m < /dev/null
AttachSpec("ShimuraQuotients.spec");
SetColumns(0);
VW := "/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/vvdata/weyl-campaign/voight-willis/";
load "/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/vvdata/weyl-campaign/voight-willis/vw3.m";
load "/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/vvdata/weyl-campaign/voight-willis/vwal.m";

D := StringToInteger(DISC); N := StringToInteger(LEVEL);
PREC := StringToInteger(PR); TARGET := StringToInteger(TGT);
BOXSZ := StringToInteger(BOX);
DN := D*N;
CC := ComplexField(PREC);

A := QuaternionAlgebra(D); O := QuaternionOrder(A, N);
G := FuchsianGroup(O);
DD := UnitDisc( : Precision := PREC);
FD := FundamentalDomain(G, DD);
rho := Max([Abs(ComplexValue(v)) : v in FD]);
H := UpperHalfPlane();
p := ComplexValue(DiscToPlane(H, DD ! 0) : Precision := PREC);
mrep := FuchsianMatrixRepresentation(A : Precision := PREC);

NT := TermsFor(rho, TARGET);
tol := ZeroLevel(rho, NT);
Q := 3*NT;

printf "=== D=%o N=%o  vol %o  rho=%o  sides %o  NT=%o  tol=%o  Q=%o  prec %o ===\n",
    D, N, RealField(5)!ArithmeticVolume(G), RealField(6)!rho, #FD, NT, RealField(3)!tol, Q, PREC;

// ---- (1) the base library on its own: Gamma relations, trivial W ----
sampG := ReducedSamples(G, O, DD, mrep, p, rho, Q, CC);
X1 := CreateShimuraQuot(D, N, {Integers() | 1});
g1, es1 := QuotientSignature(X1);
printf "Gamma: %o samples; signature (%o; %o)\n", #sampG, g1, es1;
for k in [2, 4, 6] do
    exp := DimMkSig(k, g1, es1);
    sv := SingVals(NormaliseColumns(AutSystem(sampG, k, NT, CC)));
    got, gap := GapDim(sv, tol);
    printf "  k=%o exp %o got %o%o gap %o | 6 smallest %o\n", k, exp, got,
        got eq exp select "" else "  MISMATCH", RealField(3)!gap,
        [RealField(3)!s : s in sv[1..Min(6, #sv)]];
end for;

// ---- (2) moduli of the points actually used ----
// For a Gamma sample the source is on the circle |w| = rho and the target is the reduced image.
function ModStats(samp)
    ms := Sort([Abs(s[2]) : s in samp]);
    if #ms eq 0 then return "none"; end if;
    return Sprintf("n=%o min %o med %o max %o, over rho: %o",
        #ms, RealField(4)!ms[1], RealField(4)!ms[(#ms+1) div 2], RealField(4)!ms[#ms],
        #[x : x in ms | x gt rho + 1e-10]);
end function;
printf "Gamma target moduli: %o\n", ModStats(sampG);

mu := AssociativeArray();
for m in [d : d in Divisors(DN) | d gt 1 and GCD(d, DN div d) eq 1] do
    ok, x := ALElement(O, m, BOXSZ);
    error if not ok, Sprintf("no Atkin-Lehner element of norm %o in a box of %o", m, BOXSZ);
    mu[m] := x;
    printf "  mu_%o: reduced norm %o, normalises %o\n", m, Norm(x), ALNormalises(x, O);
end for;

for m in Sort(SetToSequence(Keys(mu))) do
    sA, nbad := ALSamples(G, O, DD, mrep, p, rho, Q, CC, [mu[m]]);
    printf "AL m=%o: kept %o dropped %o; target moduli %o\n", m, #sA, nbad, ModStats(sA);
end for;

// ---- (3) does tightening a modulus threshold on the AL rows restore the gap? ----
gensAll := [d : d in Divisors(DN) | d gt 1 and GCD(d, DN div d) eq 1];
m0 := gensAll[1];
W0 := {Integers() | 1, m0};
X0 := CreateShimuraQuot(D, N, W0);
g0, es0 := QuotientSignature(X0);
sA := ALSamples(G, O, DD, mrep, p, rho, Q, CC, [mu[m0]]);
printf "\nthreshold sweep on W = {1, %o}, signature (%o; %o)\n", m0, g0, es0;
for frac in [1.0, 0.95, 0.9, 0.8, 0.7, 0.6] do
    thr := frac * rho;
    keep := [s : s in sA | Abs(s[2]) le thr and Abs(s[1]) le thr + 1e-10];
    samp := sampG cat keep;
    line := Sprintf("  thr %o*rho (%o): AL rows %o of %o |", RealField(3)!frac, RealField(4)!thr, #keep, #sA);
    for k in [2, 4] do
        exp := DimMkSig(k, g0, es0);
        sv := SingVals(NormaliseColumns(AutSystem(samp, k, NT, CC)));
        got, gap := GapDim(sv, tol);
        line cat:= Sprintf(" k=%o exp %o got %o gap %o;", k, exp, got, RealField(3)!gap);
    end for;
    print line;
end for;
printf "ALDIAG DONE\n";
quit;
