// One negative control per process.  control.m ran all four in one Magma session and the last
// two returned bit-identical residuals, because the associative arrays inside the setup record
// are shared between a record and its copy, so spoiling one spoiled the other.  Each control
// here gets its own process and its own setup, so no two can contaminate each other.
//
//   magma -b DISC:=6 K:=4 PR:=30 CUT:=40 TGT:=20 NPTS:=8 MODE:=0 control2.m < /dev/null
//
// MODE 0  correct input
//      1  one Hecke eigenvalue a_f(p) negated, at the smallest p not dividing DN
//      2  the divisor sum truncated to d = 1
//      3  R^(d) replaced by R for every d, the divisor sum otherwise intact
//      4  the eigenform of a DIFFERENT base substituted (level 10 weight 4 at D = 6)
//      5  the divisor sum missing only its TOP term d = DN, the dual lattice
//
// MODE 3 turns out NOT to be an independent control: at a ramified prime a_f(dt) = a_f(d) a_f(t),
// so with every lattice equal to R the divisor sum collapses to a scalar times the d = 1 term and
// mode 3 returns exactly TWICE mode 2, to every digit.  MODE 5 spoils the divisor sum without
// collapsing it, by dropping only the term of the dual lattice; permuting the lattices instead
// is not usable as a control because d nrd(alpha) then fails to be an integer and the run stops,
// which is itself a check that the R^(d) are the lattices the formula wants.
AttachSpec("/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/ShimuraQuotients.spec");
SetColumns(0);
load "/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/vvdata/weyl-campaign/voight-willis/vw3.m";
load "/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/vvdata/weyl-campaign/voight-willis/vwal.m";
load "/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/vvdata/weyl-campaign/nelson/nelson2.m";

D := StringToInteger(DISC); N := 1; k := StringToInteger(K);
PREC := StringToInteger(PR); cutv := StringToInteger(CUT);
TARGET := StringToInteger(TGT); npts := StringToInteger(NPTS);
mode := StringToInteger(MODE);
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

ws := [ CC | 0.40*rho*Exp(2*Pi(RR)*Ic*j/npts) + 0.09*rho*Exp(Ic*(RR!0.7)) : j in [0..npts-1] ];
zs := [ ComplexValue(DiscToPlane(H, DD ! w) : Precision := PREC) : w in ws ];
ys := [ Imaginary(z) : z in zs ];
VWvals := [ [ EvalForm(b, k, w) : w in ws ] : b in ker ];

nmax := Max(60, Ceiling(D*N*(cutv/(4*Pi(RR)))^2) + 20);
lev := mode eq 4 select 10 else D*N;
f := Newforms(CuspForms(Gamma0(lev), k))[1][1];
qq := qExpansion(f, nmax+1);
afs := [ CC ! (Rationals() ! Coefficient(qq, n)) : n in [1..nmax] ];
if mode eq 1 then
    pbad := 2; while (D*N) mod pbad eq 0 do pbad := NextPrime(pbad); end while;
    for n in [1..nmax] do
        if n mod pbad eq 0 and (n div pbad) mod pbad ne 0 then afs[n] := -afs[n]; end if;
    end for;
end if;
af := func< n | n le nmax select afs[n] else 0 >;

dat := NelsonSetup(D, N, k, PREC, cutv + 2 : OO := O, MREP := mrep);
if mode eq 2 then dat`ds := [1]; end if;
if mode eq 3 then
    b1 := dat`latb[1]; q1 := dat`latQ[1];
    for d in dat`ds do dat`latb[d] := b1; dat`latQ[d] := q1; end for;
end if;
if mode eq 5 then dat`ds := [ d : d in Divisors(D*N) | d ne D*N ]; end if;

S := [ CC | ];
for j in [1..npts] do
    s := NelsonEval(dat, zs[1], zs[j], af, RR!cutv, CC);
    Append(~S, s/(ys[1]*ys[j])^(k/2));
end for;
rows := [ S ] cat VWvals;
MM := Matrix(CC, #rows, npts,
    &cat[ [ r[j]/Sqrt(&+[Abs(x)^2 : x in r]) : j in [1..npts] ] : r in rows ]);
sv := SingValsOf(MM);
printf "MODE %o at X_0^%o(1) weight %o cut %o: residual %o   next %o\n",
       mode, D, k, cutv, RR!sv[1], RR!sv[2];
printf "  first three values: %o\n", [ CC ! x : x in S[1..3] ];
printf "CONTROL2 DONE\n";
quit;
