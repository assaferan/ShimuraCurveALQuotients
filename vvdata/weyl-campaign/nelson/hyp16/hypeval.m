// Evaluate the W-fixed weight-2 eigenforms of X_0^D(N)/W by Nelson's Theorem 3.1 at a set of
// points arranged as CIRCLES, and dump the values.
//
// WHY CIRCLES.  The canonical image of a hyperelliptic curve is a rational normal curve and
// carries no moduli at all, so the values of the differentials alone cannot give the equation:
// every product of holomorphic differentials is invariant under the hyperelliptic involution,
// and the branch points are exactly the odd information.  What does give it is ONE derivative.
// With x = u2/u1 a coordinate on the canonical conic, the function y = (dx/dz)/u1(z) is odd
// under the hyperelliptic involution and satisfies y^2 = f(x) with f of degree 8 (genus 3).
//
// A derivative taken by finite differences would throw away a third of the digits.  Taken by
// Cauchy's formula on a circle of radius rho it throws away almost none: the trapezoidal rule
// there converges like (rho/R)^NP, R the distance to the nearest pole of x, while the input
// noise is amplified only by 1/rho.  The same circle also yields the Taylor coefficients, so
// the NP/2-point estimate is available from a subset of the same values as a convergence check.
//
//   magma -b DISC:=35 WS:=1 PR:=40 CUT:=46 NC:=14 NP:=16 RHO:=0.03 OUT:=ev_35.m hypeval.m
AttachSpec("/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/ShimuraQuotients.spec");
SetColumns(0);
load "/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/vvdata/weyl-campaign/nelson/nelson3.m";

D := StringToInteger(DISC);
N := (assigned LEV) select StringToInteger(LEV) else 1;
W := { StringToInteger(x) : x in Split(WS, ",") };
PREC := StringToInteger(PR); cutv := StringToInteger(CUT);
nc := StringToInteger(NC); np := StringToInteger(NP);
CC := ComplexField(PREC); RR := RealField(PREC); Ic := CC.1;
rho := RR ! (StringToInteger(RHO) / 10000);   // RHO is given in units of 1e-4
k := 2;

g := GenusShimuraCurveQuotient(D, N, W);
printf "=== X_0^%o(%o) / %o : genus %o, cut %o, %o digits, %o circles of %o, rho %o\n",
       D, N, Sort(SetToSequence(W)), g, cutv, PREC, nc, np, RR!rho;
require g ge 3 : "genus at least 3";
require IsSquarefree(D*N) and GCD(D, N) eq 1 : "squarefree D*N assumed by the Atkin-Lehner convention";

// --- the eigenforms, exactly as stage3n.m selects them ---
nmax := Max(60, Ceiling(D*N*(cutv/(4*Pi(RR)))^2) + 20);
NFs := Newforms(CuspForms(Gamma0(D*N), k));
sel := [];
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
        if inv then Append(~sel, cs); end if;
    end for;
end for;
printf "W-fixed newforms of level exactly %o: %o (genus %o)\n", D*N, #sel, g;
require #sel eq g : "the W-fixed newforms of level exactly D*N do not span the space";

// the Atkin-Lehner eigenvalue of each form at every Hall divisor, so the analysis can tell which
// differential belongs to which quotient without re-deriving the convention
halls := [ m : m in Divisors(D*N) | GCD(m, (D*N) div m) eq 1 ];
alsig := [];
for cs in sel do
    row := [];
    for m in halls do
        e := 1;
        for p in PrimeDivisors(m) do
            e *:= (D mod p eq 0) select Round(Real(cs[p])) else -Round(Real(cs[p]));
        end for;
        Append(~row, e);
    end for;
    Append(~alsig, row);
end for;

Wtab := WkTable(k, cutv + 2, RR);
dat := NelsonSetupT(D, N, k, PREC, cutv + 2, Wtab);
sp := Sqrt(RR ! ArithmeticVolume(FuchsianGroup(dat`O)));
z1 := CC ! (0.1103 + 1.0411*Ic);

// centres spread by the square root of the arithmetic volume, as stage3.m found necessary:
// clustered points make the canonical fit report relations that are not there.
ctr := [ CC | (sp*0.73*Cos(RR ! (1.1*j^2 + 0.3*j)))
              + (0.55 + sp*0.42*(1 + Sin(RR ! (0.7*j^2 + 1.9*j))))*Ic : j in [1..nc] ];
pts := [ CC | ]; owner := []; isctr := [];
for c in [1..nc] do
    Append(~pts, ctr[c]); Append(~owner, c); Append(~isctr, 1);
    for t in [0..np-1] do
        Append(~pts, ctr[c] + rho*Exp(2*Pi(CC)*Ic*t/np));
        Append(~owner, c); Append(~isctr, 0);
    end for;
end for;
printf "%o evaluations in total\n", #pts;

F := Open(OUT, "w");
fprintf F, "// Nelson weight-2 values for X_0^%o(%o)/%o, cut %o, prec %o, rho %o\n",
        D, N, Sort(SetToSequence(W)), cutv, PREC, RR!rho;
fprintf F, "DD := %o; NN := %o; WW := %o; GG := %o; PRECV := %o; CUTV := %o;\n",
        D, N, Sort(SetToSequence(W)), g, PREC, cutv;
fprintf F, "RHO := %o; NCIRC := %o; NPC := %o;\n", RR!rho, nc, np;
fprintf F, "HALLS := %o;\n", halls;
fprintf F, "ALSIG := %o;\n", alsig;
fprintf F, "PTS := [\n";
for j in [1..#pts] do
    fprintf F, "  [%o, %o]%o\n", Real(pts[j]), Imaginary(pts[j]), j lt #pts select "," else "";
end for;
fprintf F, "];\nOWNER := %o;\nISCTR := %o;\n", owner, isctr;
fprintf F, "VALS := [\n";
delete F;

t0 := Cputime(); totpts := 0;
for j in [1..#pts] do
    S, nt := NelsonEvalMulti(dat, z1, pts[j], sel, RR!cutv, CC);
    totpts +:= nt;
    F := Open(OUT, "a");
    fprintf F, "  [%o]%o\n",
        &cat[ Sprintf("[%o, %o]%o", Real(S[i]), Imaginary(S[i]), i lt g select ", " else "")
              : i in [1..g] ],
        j lt #pts select "," else "";
    delete F;
    if j mod 10 eq 0 then printf "  %o/%o  %o s\n", j, #pts, RR!Cputime(t0); end if;
end for;
F := Open(OUT, "a");
fprintf F, "];\n";
delete F;
printf "%o evaluations: %o s, %o lattice points\n", #pts, RR!Cputime(t0), totpts;
printf "HYPEVAL DONE\n";
quit;
