// hypeval.m split across processes.  Identical point set, eigenform selection and output
// format to hypeval.m; this copy evaluates only the points JS..JE and writes
//   OUT.hdr   (the header, written only by the chunk with JS = 1)
//   OUT.pJS   (the VALS rows for its range)
// and  cat OUT.hdr OUT.p* ; echo "];"  > OUT   reassembles exactly what hypeval.m writes,
// provided the chunks tile 1..#pts.  (Rows carry their own trailing comma logic, so order of
// the parts in the cat must be by increasing JS -- use merge_parts.sh.)
//
//   magma -b DISC:=55 LEV:=7 WS:=1,35,55,77 PR:=40 CUT:=40 NC:=12 NP:=16 RHO:=300 JS:=1 JE:=51 \
//            OUT:=ev55_7_c40.m hypeval_part.m < /dev/null
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
js := StringToInteger(JS); je := StringToInteger(JE);

g := GenusShimuraCurveQuotient(D, N, W);
printf "=== X_0^%o(%o) / %o : genus %o, cut %o, %o digits, %o circles of %o, rho %o, points %o..%o\n",
       D, N, Sort(SetToSequence(W)), g, cutv, PREC, nc, np, RR!rho, js, je;
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
// The Eichler order must be the SAME in every chunk (QuaternionOrder(A, N) is not deterministic
// across processes), so it is rebuilt from the Z-basis that hypsetup.m wrote.
require assigned SETUP : "SETUP (the file hypsetup.m wrote) is required for a split evaluation";
sdat := eval (Read(SETUP) cat "\nreturn <SA, SB, OBAS>;");
Aq := QuaternionAlgebra(D);
sa, sb := StandardForm(Aq);
require sa eq sdat[1] and sb eq sdat[2] : "QuaternionAlgebra(D) came out in a different standard form than the setup file";
Oq := QuaternionOrder([ Aq ! v : v in sdat[3] ]);
require Discriminant(Oq) eq D*N : "the rebuilt order does not have discriminant D*N";
require { Eltseq(Aq ! x) : x in ZBasis(Oq) } eq { Universe([Eltseq(Aq!Aq.1)]) | v : v in sdat[3] }
     or &and[ IsCoercible(Oq, Aq ! v) : v in sdat[3] ] : "the rebuilt order is not the one in the setup file";
dat := NelsonSetupT(D, N, k, PREC, cutv + 2, Wtab : OO := Oq);
printf "Eichler order rebuilt from %o: Z-basis %o\n", SETUP, [ Eltseq(Aq ! x) : x in ZBasis(dat`O) ];
sp := Sqrt(RR ! ArithmeticVolume(FuchsianGroup(dat`O)));
z1 := CC ! (0.1103 + 1.0411*Ic);

if assigned CTRF then
    // explicit centres, one "re im" pair per line (chosen by xscan.m where the formula centres
    // leave the large-|x| end of the quartic undersampled); NC must equal the number of lines
    ctr := [ CC | ];
    for ln in Split(Read(CTRF), "\n") do
        if #ln eq 0 or ln[1] eq "/" then continue; end if;
        pr := Split(ln, " ");
        Append(~ctr, (RR ! eval(pr[1])) + (RR ! eval(pr[2]))*Ic);   // decimal literals
    end for;
    require #ctr eq nc : "NC must equal the number of centres in CTRF";
    printf "centres from %o\n", CTRF;
else
ctr := [ CC | (sp*0.73*Cos(RR ! (1.1*j^2 + 0.3*j)))
              + (0.55 + sp*0.42*(1 + Sin(RR ! (0.7*j^2 + 1.9*j))))*Ic : j in [1..nc] ];
end if;
pts := [ CC | ]; owner := []; isctr := [];
for c in [1..nc] do
    Append(~pts, ctr[c]); Append(~owner, c); Append(~isctr, 1);
    for t in [0..np-1] do
        Append(~pts, ctr[c] + rho*Exp(2*Pi(CC)*Ic*t/np));
        Append(~owner, c); Append(~isctr, 0);
    end for;
end for;
printf "%o evaluations in total, this chunk %o..%o\n", #pts, js, je;
require js ge 1 and je le #pts and js le je : "chunk range out of bounds";

if js eq 1 then
    F := Open(OUT cat ".hdr", "w");
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
end if;

pfn := OUT cat Sprintf(".p%o", js);
F := Open(pfn, "w"); delete F;
t0 := Cputime(); totpts := 0;
for j in [js..je] do
    S, nt := NelsonEvalMulti(dat, z1, pts[j], sel, RR!cutv, CC);
    totpts +:= nt;
    F := Open(pfn, "a");
    fprintf F, "  [%o]%o\n",
        &cat[ Sprintf("[%o, %o]%o", Real(S[i]), Imaginary(S[i]), i lt g select ", " else "")
              : i in [1..g] ],
        j lt #pts select "," else "";
    delete F;
    if (j - js + 1) mod 5 eq 0 then printf "  %o/%o  %o s\n", j - js + 1, je - js + 1, RR!Cputime(t0); end if;
end for;
printf "%o evaluations: %o s, %o lattice points\n", je - js + 1, RR!Cputime(t0), totpts;
printf "HYPEVAL PART DONE\n";
quit;
