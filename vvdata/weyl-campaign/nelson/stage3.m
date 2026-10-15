// Stage 3.  Max Noether's count of quadrics through the canonical image, with the basis of
// differentials supplied by Nelson's evaluator instead of by q-expansions or a function field.
//
// The differentials of X_0^D(N)/W are the weight-2 Hecke eigenforms fixed by W, and Theorem 3.1
// evaluates each of them from its eigenvalues.  Fixing one point z_1 and reading
// conj(f_i(z_1)) f_i(z_j) gives each form with its own nonzero constant, which is a change of
// coordinate on the canonical image and changes no rank.  A genus-g curve has g(g+1)/2 quadratic
// monomials in those coordinates, and the number of relations among their values is
// (g-2)(g-3)/2 if the curve is not hyperelliptic and (g-1)(g-2)/2 if it is -- 0 against 1 at
// genus 3.
//
//   magma -b DISC:=35 WS:=1 PR:=30 CUT:=36 NPTS:=14 stage3.m < /dev/null
AttachSpec("/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/ShimuraQuotients.spec");
SetColumns(0);
load "/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/vvdata/weyl-campaign/nelson/nelson3.m";

D := StringToInteger(DISC); N := 1;
W := { StringToInteger(x) : x in Split(WS, ",") };
PREC := StringToInteger(PR); cutv := StringToInteger(CUT); npts := StringToInteger(NPTS);
CC := ComplexField(PREC); RR := RealField(PREC); Ic := CC.1;
k := 2;

g := GenusShimuraCurveQuotient(D, N, W);
printf "=== X_0^%o(%o) / %o : genus %o, cut %o, %o points, %o digits\n",
       D, N, Sort(SetToSequence(W)), g, cutv, npts, PREC;
require g ge 3 : "Noether's count needs genus at least 3";

// The eigenforms.  At N = 1 and D squarefree every form of the D-new space of level D is a
// newform of level exactly D, which is what Theorem 3.1 is stated for.  On the quaternion side
// the Atkin-Lehner eigenvalue at p | D is MINUS the classical one, so the eigenvalue of w_m is
// the product of the classical a_f(p) over p | m; the count against the genus below is what
// checks that convention rather than assuming it.
nmax := Max(60, Ceiling(D*N*(cutv/(4*Pi(RR)))^2) + 20);
NFs := Newforms(CuspForms(Gamma0(D*N), k));
sel := []; lab := [];
for orb in NFs do
    for f in orb do
        KF := BaseRing(Parent(f));
        qq := qExpansion(f, nmax+1);
        if Degree(KF) eq 1 then
            cs := [ CC ! (Rationals() ! Coefficient(qq, n)) : n in [1..nmax] ];
        else
            cj := Conjugates(KF.1 : Precision := PREC+10);
            // one embedding per conjugate: pick the one matching this member of the orbit
            idx := Index(orb, f);
            cs := [ CC | Evaluate(Polynomial(Rationals(), Eltseq(Coefficient(qq,n))), cj[idx])
                    : n in [1..nmax] ];
        end if;
        inv := true;
        for m in W do
            if m eq 1 then continue; end if;
            e := &*[ cs[p] : p in PrimeDivisors(m) ];
            if Abs(e - 1) gt 0.01 then inv := false; end if;
        end for;
        if inv then
            Append(~sel, cs);
            Append(~lab, [ RR ! Real(cs[n]) : n in [2..6] ]);
        end if;
    end for;
end for;
printf "eigenforms fixed by W: %o (genus says %o)\n", #sel, g;
require #sel eq g : "the number of W-fixed eigenforms does not match the genus";
for l in lab do printf "   a_2..a_6 = %o\n", l; end for;

// Evaluate all g forms at npts points, in ONE lattice enumeration per point.
Wtab := WkTable(k, cutv + 2, RR);
dat := NelsonSetupT(D, N, k, PREC, cutv + 2, Wtab);
// The points must be spread widely enough that no two are close modulo the group: at a base
// with a large fundamental domain a cluster of nearby points gives rows that are nearly
// proportional, and the monomial matrix is then rank deficient for a reason that has nothing to
// do with the canonical image.  The spread is scaled by the arithmetic volume.
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

// Noether's matrix: one row per point, one column per quadratic monomial.
mons := [ [a,b] : b in [a..g], a in [1..g] ];
M := Matrix(CC, npts, #mons,
    &cat[ [ vals[j][m[1]]*vals[j][m[2]] : m in mons ] : j in [1..npts] ]);
// equilibrate rows and columns so the singular values mean something
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

// the number of relations is the number of singular values at the noise level
best := 1; bestr := RR ! 1;
for i in [1..#sv-1] do
    if sv[i] gt 0 and sv[i+1]/sv[i] gt bestr then bestr := sv[i+1]/sv[i]; best := i; end if;
end for;
// A verdict is only allowed on a CLEAN gap.  Reading off the largest ratio always names a
// number, even on a smooth ramp where the matrix is simply ill conditioned, which is how a
// spurious relation gets reported.  MINGAP is the ratio the gap has to reach.
mingap := RR ! 10^4;
if sv[1] gt 10^(-4) then
    nrel := 0; ok := true;              // nothing near zero at all: no relations
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
printf "STAGE3 DONE\n";
quit;
