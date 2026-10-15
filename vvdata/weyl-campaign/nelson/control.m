// Negative controls for the stage-1 agreement.
//
// The agreement between Nelson's evaluator and the Voight-Willis expansions is only evidence if
// it FAILS when the input is wrong, and fails for the stated reason: the eigenvalues are what
// picks out the form, so changing one of them must stop the values lying in the Voight-Willis
// space.  Three inputs are spoiled in turn, one at a time:
//
//   1. one Hecke eigenvalue a_f(p) at a prime p not dividing DN;
//   2. the coefficients c_f(d) of the divisor sum, by dropping every d > 1;
//   3. the lattices, by using R in place of R^(d) throughout.
//
//   magma -b DISC:=6 K:=4 PR:=30 CUT:=40 TGT:=20 NPTS:=8 control.m < /dev/null
AttachSpec("/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/ShimuraQuotients.spec");
SetColumns(0);
load "/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/vvdata/weyl-campaign/voight-willis/vw3.m";
load "/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/vvdata/weyl-campaign/voight-willis/vwal.m";
load "/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/vvdata/weyl-campaign/nelson/nelson2.m";

D := StringToInteger(DISC); N := 1; k := StringToInteger(K);
PREC := StringToInteger(PR); cutv := StringToInteger(CUT);
TARGET := StringToInteger(TGT); npts := StringToInteger(NPTS);
CC := ComplexField(PREC); RR := RealField(PREC); Ic := CC.1;

function KernelVecs(Mx, r)
    n := Ncols(Mx);
    cn := [ Sqrt(&+[ Abs(Mx[i][j])^2 : i in [1..Nrows(Mx)] ]) : j in [1..n] ];
    Sig, L, Rt := SingularValueDecomposition(RealForm(NormaliseColumns(Mx)));
    rr := Min(Nrows(Sig), Ncols(Sig));
    return [ [ (Rt[rr-2*(s-1)][i] + Ic*Rt[rr-2*(s-1)][n+i])/cn[i] : i in [1..n] ] : s in [1..r] ];
end function;
function SingValsOf(Mx)
    MR := RealForm(Mx); Sig := SingularValueDecomposition(MR);
    r := Min(Nrows(MR), Ncols(MR));
    sv := Reverse([ Abs(Sig[i][i]) : i in [1..r] ]);
    return [ sv[i] : i in [1..#sv by 2] ];
end function;

A := QuaternionAlgebra(D); O := QuaternionOrder(A, N);
G := FuchsianGroup(O); H := UpperHalfPlane();
_, DD, FD := RecentredDisc(G, H, PREC, 3);
p0 := ComplexValue(DiscToPlane(H, DD ! 0) : Precision := PREC);
mrep := FuchsianMatrixRepresentation(A : Precision := PREC);
rho := Max([Abs(ComplexValue(v)) : v in FD]);
NT := TermsFor(rho, TARGET);
samp := ReducedSamples(G, O, DD, mrep, p0, rho, 3*NT, CC);
Mx := AutSystem(samp, k, NT, CC);
X1 := CreateShimuraQuot(D, N, {Integers() | 1});
g1, es1 := QuotientSignature(X1);
dimexp := DimMkSig(k, g1, es1);
ker := KernelVecs(Mx, dimexp);
printf "=== X_0^%o(1) weight %o, dim %o, rho %o, NT %o, cut %o\n", D, k, dimexp, RR!rho, NT, cutv;

ws := [ CC | 0.40*rho*Exp(2*Pi(RR)*Ic*j/npts) + 0.09*rho*Exp(Ic*(RR!0.7)) : j in [0..npts-1] ];
zs := [ ComplexValue(DiscToPlane(H, DD ! w) : Precision := PREC) : w in ws ];
ys := [ Imaginary(z) : z in zs ];
VWvals := [ [ EvalForm(b, k, w) : w in ws ] : b in ker ];

f := Newforms(CuspForms(Gamma0(D*N), k))[1][1];
nmax := Max(60, Ceiling(D*N*(cutv/(4*Pi(RR)))^2) + 20);
qq := qExpansion(f, nmax+1);
afs := [ CC ! (Rationals() ! Coefficient(qq, n)) : n in [1..nmax] ];

// the smallest prime not dividing DN, and the coefficients with a_f at that prime negated
pbad := 2; while (D*N) mod pbad eq 0 do pbad := NextPrime(pbad); end while;
bad := afs;
for n in [1..nmax] do
    if n mod pbad eq 0 and (n div pbad) mod pbad ne 0 then bad[n] := -bad[n]; end if;
end for;

function Residual(dat, af, ds)
    S := [ CC | ];
    for j in [1..npts] do
        dd := dat; dd`ds := ds;
        s := NelsonEval(dd, zs[1], zs[j], af, RR!cutv, CC);
        Append(~S, s/(ys[1]*ys[j])^(k/2));
    end for;
    rows := [ S ] cat VWvals;
    MM := Matrix(CC, #rows, npts,
        &cat[ [ r[j]/Sqrt(&+[Abs(x)^2 : x in r]) : j in [1..npts] ] : r in rows ]);
    return SingValsOf(MM);
end function;

dat := NelsonSetup(D, N, k, PREC, cutv + 2 : OO := O, MREP := mrep);
afgood := func< n | n le nmax select afs[n] else 0 >;
afbad  := func< n | n le nmax select bad[n] else 0 >;

s := Residual(dat, afgood, dat`ds);
printf "  correct input                          : residual %o   (next %o)\n", RR!s[1], RR!s[2];
s := Residual(dat, afbad, dat`ds);
printf "  a_f(%o) negated                         : residual %o   (next %o)\n", pbad, RR!s[1], RR!s[2];
s := Residual(dat, afgood, [1]);
printf "  divisor sum truncated to d = 1 only    : residual %o   (next %o)\n", RR!s[1], RR!s[2];
dflat := dat;
for d in dat`ds do dflat`latb[d] := dat`latb[1]; dflat`latQ[d] := dat`latQ[1]; end for;
s := Residual(dflat, afgood, dat`ds);
printf "  R^(d) replaced by R for every d        : residual %o   (next %o)\n", RR!s[1], RR!s[2];
printf "CONTROL DONE\n";
quit;
