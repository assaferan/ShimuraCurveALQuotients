// Stage 3 at general N.  Identical to stage3.m except that the level is not forced to 1:
//
//  * only the newforms of level EXACTLY D*N are used, which is what Theorem 3.1 is stated for;
//    nold.m is what says whether those span the W-fixed space, and this script refuses to run
//    unless their number equals the genus, so a base short of the space stops here rather than
//    returning a count taken over a subspace;
//  * the Atkin-Lehner eigenvalue on the quaternion side is minus the classical one at p | D and
//    equal to it at p | N, and the classical eigenvalue of w_p in weight 2 at squarefree level is
//    -a_f(p); so the quaternionic eigenvalue of w_m is the product over p | m of a_f(p) at p | D
//    and of -a_f(p) at p | N.  At N = 1 this is the convention stage3.m uses.
//
//   magma -b DISC:=55 LEV:=7 WS:=1,35,55,77 PR:=30 CUT:=40 NPTS:=18 stage3n.m < /dev/null
AttachSpec("/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/ShimuraQuotients.spec");
SetColumns(0);
load "/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/vvdata/weyl-campaign/nelson/nelson3hi.m";

D := StringToInteger(DISC);
N := (assigned LEV) select StringToInteger(LEV) else 1;
W := { StringToInteger(x) : x in Split(WS, ",") };
PREC := StringToInteger(PR); cutv := StringToInteger(CUT); npts := StringToInteger(NPTS);
CC := ComplexField(PREC); RR := RealField(PREC); Ic := CC.1;
k := 2;

g := GenusShimuraCurveQuotient(D, N, W);
printf "=== X_0^%o(%o) / %o : genus %o, cut %o, %o points, %o digits\n",
       D, N, Sort(SetToSequence(W)), g, cutv, npts, PREC;
require g ge 3 : "Noether's count needs genus at least 3";
require IsSquarefree(D*N) and GCD(D, N) eq 1 : "the Atkin-Lehner convention below assumes D*N squarefree";

nmax := Max(60, Ceiling(D*N*(cutv/(4*Pi(RR)))^2) + 20);
NFs := Newforms(CuspForms(Gamma0(D*N), k));
sel := []; lab := [];
for orb in NFs do
    for f in orb do
        if Level(Parent(f)) ne D*N then continue; end if;
        KF := BaseRing(Parent(f));
        qq := qExpansion(f, nmax+1);
        if Degree(KF) eq 1 then
            cs := [ CC ! (Rationals() ! Coefficient(qq, n)) : n in [1..nmax] ];
        else
            cj := Conjugates(KF.1 : Precision := PREC+10);
            idx := Index(orb, f);
            cs := [ CC | Evaluate(Polynomial(Rationals(), Eltseq(Coefficient(qq,n))), cj[idx])
                    : n in [1..nmax] ];
        end if;
        inv := true;
        for m in W do
            if m eq 1 then continue; end if;
            e := CC ! 1;
            for p in PrimeDivisors(m) do
                e *:= (D mod p eq 0) select cs[p] else -cs[p];
            end for;
            if Abs(e - 1) gt 0.01 then inv := false; end if;
        end for;
        if inv then
            Append(~sel, cs);
            Append(~lab, [ RR ! Real(cs[n]) : n in [2..6] ]);
        end if;
    end for;
end for;
printf "eigenforms fixed by W, of level exactly %o: %o (genus says %o)\n", D*N, #sel, g;
require #sel eq g :
    "the W-fixed newforms of level exactly D*N do not span the space (see nold.m); Theorem 3.1 does not reach this base";
for l in lab do printf "   a_2..a_6 = %o\n", l; end for;

Wtab := WkTable(k, cutv + 2, RR);
dat := NelsonSetupT(D, N, k, PREC, cutv + 2, Wtab);
sp := Sqrt(RR ! ArithmeticVolume(FuchsianGroup(dat`O)));
z1 := CC ! (0.1103 + 1.0411*Ic);
zs := [ z1 ] cat [ CC | (sp*0.73*Cos(RR ! (1.1*j^2 + 0.3*j)))
                        + (0.55 + sp*0.42*(1 + Sin(RR ! (0.7*j^2 + 1.9*j))))*Ic
                   : j in [1..npts-1] ];
vals := [];
t0 := Cputime(); totpts := 0;
for j in [1..npts] do
    S, nt := NelsonEvalMulti(dat, z1, zs[j], sel, RR!cutv, CC);
    Append(~vals, S);
    totpts +:= nt;
end for;
printf "%o evaluations of %o forms: %o s, %o lattice points total\n",
       npts, g, RR!Cputime(t0), totpts;

mons := [ [a,b] : b in [a..g], a in [1..g] ];
M := Matrix(CC, npts, #mons,
    &cat[ [ vals[j][m[1]]*vals[j][m[2]] : m in mons ] : j in [1..npts] ]);
for i in [1..npts] do
    nm := Sqrt(&+[ Abs(M[i][c])^2 : c in [1..#mons] ]);
    for c in [1..#mons] do M[i][c] := M[i][c]/nm; end for;
end for;
for c in [1..#mons] do
    nm := Sqrt(&+[ Abs(M[i][c])^2 : i in [1..npts] ]);
    for i in [1..npts] do M[i][c] := M[i][c]/nm; end for;
end for;
MR := Matrix(RR, 2*npts, 2*#mons, [ RR ! 0 : i in [1..4*npts*#mons] ]);
for i in [1..npts] do for j in [1..#mons] do
    a := Real(M[i][j]); b := Imaginary(M[i][j]);
    MR[i][j] := a; MR[i][#mons+j] := -b; MR[npts+i][j] := b; MR[npts+i][#mons+j] := a;
end for; end for;
Sig := SingularValueDecomposition(MR);
r := Min(Nrows(MR), Ncols(MR));
sv := Reverse([ Abs(Sig[i][i]) : i in [1..r] ]);
sv := [ sv[i] : i in [1..#sv by 2] ];
printf "singular values of the %o x %o monomial matrix:\n  %o\n", npts, #mons,
       [ RR ! x : x in sv ];
printf "consecutive ratios: %o\n", [ RR ! (sv[i+1]/sv[i]) : i in [1..#sv-1] ];

best := 1; bestr := RR ! 1;
for i in [1..#sv-1] do
    if sv[i] gt 0 and sv[i+1]/sv[i] gt bestr then bestr := sv[i+1]/sv[i]; best := i; end if;
end for;
mingap := RR ! 10^4;
if sv[1] gt 10^(-4) then
    nrel := 0; ok := true;
else
    nrel := best; ok := bestr ge mingap;
end if;
printf "largest ratio jump %o after %o value(s); clean gap required %o\n", RR!bestr, best, mingap;
if not ok then
    printf "=> NO CLEAN GAP: the spectrum is a ramp, the count does not resolve\n";
    printf "VERDICT: INCONCLUSIVE\n";
else
    printf "=> quadric relations: %o   (hyperelliptic would be %o, non-hyperelliptic %o)\n",
           nrel, (g-1)*(g-2) div 2, (g-2)*(g-3) div 2;
    printf "VERDICT: %o\n",
       nrel eq (g-1)*(g-2) div 2 select "HYPERELLIPTIC"
       else (nrel eq (g-2)*(g-3) div 2 select "NOT hyperelliptic" else "INCONCLUSIVE");
end if;
printf "STAGE3N DONE\n";
quit;
