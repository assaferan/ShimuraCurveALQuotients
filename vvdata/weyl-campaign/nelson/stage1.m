// Stage 1.  Validate Nelson's evaluator against the Voight-Willis expansions at the same points.
//
// For each Hecke eigenform of weight k on X_0^D(1) -- specified by the eigenvalues of the
// corresponding newform of level D, which is the only input Nelson's formula takes -- evaluate it
// at m points by Theorem 3.1, and ask whether the resulting vector of values lies in the span of
// the Voight-Willis basis of the same space evaluated at the same points.  With m points and a
// space of dimension g that is m - g independent predictions per eigenform.
//
//   magma -b DISC:=6 K:=4 PR:=30 CUT:=50 TGT:=20 NPTS:=8 stage1.m < /dev/null
AttachSpec("/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/ShimuraQuotients.spec");
SetColumns(0);
load "/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/vvdata/weyl-campaign/voight-willis/vw3.m";
load "/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/vvdata/weyl-campaign/voight-willis/vwal.m";
load "/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/vvdata/weyl-campaign/nelson/nelson2.m";

D := StringToInteger(DISC); N := 1; k := StringToInteger(K);
PREC := StringToInteger(PR); cutv := StringToInteger(CUT);
TARGET := StringToInteger(TGT); npts := StringToInteger(NPTS);
CC := ComplexField(PREC); RR := RealField(PREC); Ic := CC.1;

// --- the smallest r right singular vectors, undoing the column equilibration
function KernelVecs(Mx, r)
    n := Ncols(Mx);
    cn := [ Sqrt(&+[ Abs(Mx[i][j])^2 : i in [1..Nrows(Mx)] ]) : j in [1..n] ];
    MN := NormaliseColumns(Mx);
    MR := RealForm(MN);
    Sig, L, Rt := SingularValueDecomposition(MR);
    rr := Min(Nrows(MR), Ncols(MR));
    out := [];
    for s in [1..r] do
        j := rr - 2*(s-1);                                 // each complex value appears twice
        v := [ CC | Rt[j][i] + Ic*Rt[j][n+i] : i in [1..n] ];
        Append(~out, [ v[i]/cn[i] : i in [1..n] ]);
    end for;
    return out;
end function;

function SingValsOf(Mx)
    Sig := SingularValueDecomposition(RealForm(Mx));
    r := Min(Nrows(RealForm(Mx)), Ncols(RealForm(Mx)));
    sv := Reverse([ Abs(Sig[i][i]) : i in [1..r] ]);
    return [ sv[i] : i in [1..#sv by 2] ];
end function;

A := QuaternionAlgebra(D); O := QuaternionOrder(A, N);
G := FuchsianGroup(O);
H := UpperHalfPlane();
_, DD, FD := RecentredDisc(G, H, PREC, 3);
p := ComplexValue(DiscToPlane(H, DD ! 0) : Precision := PREC);
mrep := FuchsianMatrixRepresentation(A : Precision := PREC);
rho := Max([Abs(ComplexValue(v)) : v in FD]);
NT := TermsFor(rho, TARGET);
printf "=== X_0^%o(1), weight %o:  rho=%o  NT=%o  prec=%o  cut=%o\n",
       D, k, RR!rho, NT, PREC, cutv;

// ---- Voight-Willis: a basis of the space
t0 := Cputime();
samp := ReducedSamples(G, O, DD, mrep, p, rho, 3*NT, CC);
Mx := AutSystem(samp, k, NT, CC);
sv := SingVals(NormaliseColumns(Mx));
X1 := CreateShimuraQuot(D, N, {Integers() | 1});
g1, es1 := QuotientSignature(X1);
dimexp := DimMkSig(k, g1, es1);
printf "VW: %o samples, dim expected %o, smallest svs %o  (%o s)\n",
       #samp, dimexp, [RR!x : x in sv[1..Min(5,#sv)]], RR!Cputime(t0);
ker := KernelVecs(Mx, dimexp);

// ---- the sample points
ws := [ CC | 0.40*rho*Exp(2*Pi(RR)*Ic*j/npts) + 0.09*rho*Exp(Ic*(RR!0.7)) : j in [0..npts-1] ];
zs := [ ComplexValue(DiscToPlane(H, DD ! w) : Precision := PREC) : w in ws ];
ys := [ Imaginary(z) : z in zs ];
VWvals := [ [ EvalForm(b, k, w) : w in ws ] : b in ker ];

// ---- Nelson: one evaluation per eigenform
NFs := Newforms(CuspForms(Gamma0(D*N), k));
printf "newform orbits of level %o weight %o: %o, degrees %o\n",
       D*N, k, #NFs, [ Degree(BaseRing(Parent(orb[1]))) : orb in NFs ];
dat := NelsonSetup(D, N, k, PREC, cutv + 2 : OO := O, MREP := mrep);

nform := 0;
for orb in NFs do
  for f in orb do                             // every conjugate is its own eigenform
    nform +:= 1;
    nmax := Max(60, Ceiling(D*N*(cutv/(4*Pi(RR)))^2) + 20);
    qq := qExpansion(f, nmax+1);
    KF := BaseRing(Parent(f));
    if Degree(KF) eq 1 then
        afs := [ CC ! (Rationals() ! Coefficient(qq, n)) : n in [1..nmax] ];
    else
        cj := Conjugates(KF.1 : Precision := PREC)[1];
        afs := [ CC | Evaluate(Polynomial(Eltseq(Coefficient(qq,n))), cj) : n in [1..nmax] ];
    end if;
    af := func< n | n le nmax select afs[n] else 0 >;
    t1 := Cputime();
    S := [ CC | ]; nt := 0; nmx := 0;
    for j in [1..npts] do
        s, c, m := NelsonEval(dat, zs[1], zs[j], af, RR!cutv, CC);
        Append(~S, s / (ys[1]*ys[j])^(k/2));              // = conj(f(z1)) f(z_j)
        nt +:= c; nmx := Max(nmx, m);
    end for;
    tN := Cputime(t1);
    rows := [ S ] cat VWvals;
    MM := Matrix(CC, #rows, npts, &cat[ [ r[j]/Sqrt(&+[Abs(x)^2 : x in r]) : j in [1..npts] ] : r in rows ]);
    s2 := SingValsOf(MM);
    printf "\n form %o (a_2..a_7 = %o):\n", nform, [ RR ! Real(af(n)) : n in [2..7] ];
    printf "   Nelson: %o s for %o points, %o lattice points, a_f needed to n=%o\n",
           RR!tN, npts, nt, nmx;
    printf "   values conj(f(z1))f(z_j): %o\n", [ CC ! x : x in S[1..Min(3,npts)] ];
    printf "   rank test on [Nelson ; VW basis] (%o x %o): singular values %o\n",
           #rows, npts, [RR!x : x in s2];
    printf "   => residual %o, next %o, GAP %o  (%o independent predictions)\n",
           RR!s2[1], RR!s2[2], RR!(s2[2]/s2[1]), npts - dimexp;
  end for;
end for;
printf "\nSTAGE1 DONE\n";
quit;
