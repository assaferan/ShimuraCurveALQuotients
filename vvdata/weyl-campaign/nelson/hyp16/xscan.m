// Cheap scan of |x| over a grid of points, to choose extra circle centres where |x| is large.
//
// WHY.  At X_0^6(185) the formula centres of hypeval.m gave |x| in 0.1..0.25 (relative units)
// at ten of twelve centres, so the quartic in t = x^2 was fitted on data spanning two decades of
// t with almost all of it at the small end; the coefficients of t^3, t^4 rest on two points and
// the invariants B, C come out to 1e-6 while A comes out to 1e-11.  The remedy is more centres
// at large |x|.  x = w2/w1 is read off the canonical conic of an existing dump (same Eichler
// order, same z1, same eigenform normalisation), so |x| here is on the SAME scale as that dump.
//
//   magma -b IN:=ev_6_185_c40.m SETUP:=setup_6_185.m CUT:=24 PR:=30 GX:=10 GY:=8 OUT:=xscan_6_185.txt xscan.m
AttachSpec("/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/ShimuraQuotients.spec");
SetColumns(0);
load "/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/vvdata/weyl-campaign/nelson/nelson3.m";
dmp := eval (Read(IN) cat
   "\nreturn <DD, NN, WW, GG, PRECV, CUTV, RHO, NCIRC, NPC, HALLS, ALSIG, PTS, OWNER, ISCTR, VALS>;");
D := dmp[1]; N := dmp[2]; W := Set(dmp[3]); g := dmp[4]; PRECV := dmp[5];
HALLS := dmp[10]; ALSIG := dmp[11]; PTS := dmp[12]; ISCTR := dmp[14]; VALS := dmp[15]; np := dmp[9];
CCd := ComplexField(PRECV); RRd := RealField(PRECV); Icd := CCd.1; R6 := RealField(6);
pts := [ CCd | PTS[j][1] + Icd*PTS[j][2] : j in [1..#PTS] ];
Gv := [ [ CCd | (VALS[j][i][1] + Icd*VALS[j][i][2]) / Imaginary(pts[j]) : i in [1..g] ] : j in [1..#PTS] ];
// --- conic and normal form, as hypfit2.m ---
uix := 0;
for i in [1..#HALLS] do
    if HALLS[i] eq 1 then continue; end if;
    if #[ r : r in ALSIG | r[i] eq 1 ] eq 1 then uix := i; break; end if;
end for;
pind := [ j : j in [1..g] | ALSIG[j][uix] eq 1 ][1];
aind := [ j : j in [1..g] | ALSIG[j][uix] eq -1 ];
idx := [ j : j in [1..#PTS] | ISCTR[j] eq 1 ] cat [ j : j in [1..#PTS] | ISCTR[j] eq 0 and (j mod np) eq 3 ];
mons := [ [a,b] : b in [a..g], a in [1..g] ];
Mx := Matrix(CCd, #idx, #mons, &cat[ [ Gv[j][m[1]]*Gv[j][m[2]] : m in mons ] : j in idx ]);
for i in [1..#idx] do
    nm := Sqrt(&+[ Abs(Mx[i][c])^2 : c in [1..#mons] ]);
    for c in [1..#mons] do Mx[i][c] := Mx[i][c]/nm; end for;
end for;
cn := [ Sqrt(&+[ Abs(Mx[i][c])^2 : i in [1..#idx] ]) : c in [1..#mons] ];
for c in [1..#mons] do for i in [1..#idx] do Mx[i][c] := Mx[i][c]/cn[c]; end for; end for;
n2 := #idx;
MR := Matrix(RRd, 2*n2, 2*#mons, [ RRd ! 0 : i in [1..4*n2*#mons] ]);
for i in [1..n2] do for j in [1..#mons] do
    a := Real(Mx[i][j]); b := Imaginary(Mx[i][j]);
    MR[i][j] := a; MR[i][#mons+j] := -b; MR[n2+i][j] := b; MR[n2+i][#mons+j] := a;
end for; end for;
Sg, U, V := SingularValueDecomposition(MR);
kcol := Ncols(V);
function KernRes(vec)
    s := RRd ! 0;
    for i in [1..#idx] do s +:= Abs(&+[ Mx[i][c]*vec[c] : c in [1..#mons] ])^2; end for;
    return Sqrt(s)/Sqrt(&+[ Abs(c)^2 : c in vec ]);
end function;
vcol := [ CCd | V[j][kcol] + Icd*V[#mons+j][kcol] : j in [1..#mons] ];
vrow := [ CCd | V[kcol][j] + Icd*V[kcol][#mons+j] : j in [1..#mons] ];
kvec := KernRes(vcol) lt KernRes(vrow) select vcol else vrow;
error if Minimum(KernRes(vcol), KernRes(vrow)) gt 10^(-8), "no conic in the dump";
qc := [ CCd | kvec[j]/cn[j] : j in [1..#mons] ];
Q := ZeroMatrix(CCd, 3, 3);
for c in [1..#mons] do
    a := mons[c][1]; b := mons[c][2];
    if a eq b then Q[a][a] := qc[c]; else Q[a][b] := qc[c]/2; Q[b][a] := qc[c]/2; end if;
end for;
a1 := aind[1]; a2 := aind[2];
qaa := Q[a1][a1]; qab := Q[a1][a2];
r1 := (-qab + Sqrt(qab^2 - qaa*Q[a2][a2]))/qaa;
function Xof(v)   // v = the g values (already divided by Im z)
    return v[pind]/(v[a1] - r1*v[a2]);
end function;
xc := [ Abs(Xof(Gv[j])) : j in [1..#PTS] | ISCTR[j] eq 1 ];
sigx := Exp((&+[ Log(x) : x in xc ])/#xc);
printf "existing centres: |x|/sigx = %o  (sigx = %o)\n", [ R6!(x/sigx) : x in xc ], R6!sigx;

// --- the scan ---
PREC := StringToInteger(PR); cutv := StringToInteger(CUT); gx := StringToInteger(GX); gy := StringToInteger(GY);
CC := ComplexField(PREC); RR := RealField(PREC); Ic := CC.1;
k := 2;
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
            cs := [ CC | Evaluate(Polynomial(Rationals(), Eltseq(Coefficient(qq,n))), cj[Index(orb, f)]) : n in [1..nmax] ];
        end if;
        inv := true;
        for m in W do
            if m eq 1 then continue; end if;
            e := CC ! 1;
            for p in PrimeDivisors(m) do e *:= (D mod p eq 0) select cs[p] else -cs[p]; end for;
            if Abs(e - 1) gt 0.01 then inv := false; end if;
        end for;
        if inv then Append(~sel, cs); end if;
    end for;
end for;
require #sel eq g : "newforms do not span";
sdat := eval (Read(SETUP) cat "\nreturn <SA, SB, OBAS>;");
Aq := QuaternionAlgebra(D);
Oq := QuaternionOrder([ Aq ! v : v in sdat[3] ]);
require Discriminant(Oq) eq D*N : "order";
Wtab := WkTable(k, cutv + 2, RR);
dat := NelsonSetupT(D, N, k, PREC, cutv + 2, Wtab : OO := Oq);
sp := Sqrt(RR ! ArithmeticVolume(FuchsianGroup(dat`O)));
z1 := CC ! (0.1103 + 1.0411*Ic);
cands := [ CC | ];
// default rectangle = the span of hypeval.m's formula centres; X0,X1,Y0,Y1 override it
x0 := assigned X0 select RR ! eval(X0) else -0.73*sp;  x1 := assigned X1 select RR ! eval(X1) else 0.73*sp;
y0 := assigned Y0 select RR ! eval(Y0) else RR ! 0.55; y1 := assigned Y1 select RR ! eval(Y1) else 0.55 + 0.84*sp;
for i in [1..gx] do for j in [1..gy] do
    Append(~cands, (x0 + (x1 - x0)*(i - 0.5)/gx) + (y0 + (y1 - y0)*(j - 0.5)/gy)*Ic);
end for; end for;
// consistency: re-evaluate the first two existing centres at this cutoff and compare |x|
cands := [ CC | CC!pts[1], CC!pts[np+2] ] cat cands;
printf "%o candidates (first two are existing centres 1 and 2 for comparison: |x|/sigx %o, %o)\n",
       #cands, R6!(xc[1]/sigx), R6!(xc[2]/sigx);
F := Open(OUT, "w"); fprintf F, "// re im |x|/sigx   (sigx = %o)\n", sigx; delete F;
t0 := Cputime();
for j in [1..#cands] do
    S, nt := NelsonEvalMulti(dat, z1, cands[j], sel, RR!cutv, CC);
    v := [ CCd | (CCd!S[i]) / Imaginary(CCd!cands[j]) : i in [1..g] ];
    ax := Abs(Xof(v))/sigx;
    F := Open(OUT, "a");
    fprintf F, "%o %o %o\n", Real(cands[j]), Imaginary(cands[j]), R6!ax;
    delete F;
    if j le 2 or j mod 10 eq 0 then printf "  %o: z = %o, |x|/sigx = %o   (%o s)\n", j, cands[j], R6!ax, R6!Cputime(t0); end if;
end for;
print "XSCAN DONE";
quit;
