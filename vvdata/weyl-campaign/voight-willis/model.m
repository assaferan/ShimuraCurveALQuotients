// Build the double cover of the star curve for one Atkin-Lehner subgroup, from the expansions.
//
//   magma -b DISC:=35 LEVEL:=1 WSET:=1,5 PR:=80 TGT:=20 QMUL:=3 BOX:=3 model.m < /dev/null
//
// THE SHAPE OF THE ANSWER.  For W of index 2 in the full Atkin-Lehner group W*, the quotient
// X/W is a double cover of the star curve X* = X/W*, so when the star has genus 0 the model is
// y^2 = f(t) over a hauptmodul t of the star.  Two functions have to come out of the expansions:
//   t = f1/f2,  a ratio of two forms on the star spanning a pencil with an imposed base locus;
//   y = g/h,    g changing sign under the deck transformation, h a form on the star, same weight.
// y^2 is invariant under all of W*, hence a rational function of t, and fitting P(t) - y^2 Q(t)
// recovers it.  The branch points are the roots of P*Q of odd multiplicity, and they are what can
// be compared with a published equation, since the two coordinates on the star differ by a
// Mobius change and the branch sets must correspond under one.
//
// WHERE THE FORMS ARE EXPANDED.  On the domain of the NORMALISER, not of Gamma.  The term count
// is NT = T ln10 / (-ln rho) and the decomposition costs NT^3, so the circumradius of the domain
// the points reduce into is the only number that matters for cost; the normaliser's domain is
// 1/#W* of the area and has a smaller one.  At X_0^10(9) that is 0.837 against 0.955, which is
// 259 terms against 1004 and a factor of forty-five in the decomposition.  Each sample is
// labelled by the coset of the normaliser its reduction landed in, so one set of samples cuts out
// the invariant, star and anti-invariant spaces at once.
AttachSpec("/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/ShimuraQuotients.spec");
SetColumns(0);
load "/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/vvdata/weyl-campaign/voight-willis/vw3.m";
load "/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/vvdata/weyl-campaign/voight-willis/vwal.m";
load "/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/vvdata/weyl-campaign/voight-willis/vwmodel.m";

// Magma buffers stdout when it is a file, so a log that has not grown says nothing about whether
// a run is alive.  Every stage also appends one line to a progress file, closed each time and so
// flushed.  ⚠ Magma reads a script as it runs it, so this file must not be edited while a run of
// it is in flight -- doing so produced a syntax error in the middle of a finished computation.
PROG := assigned PROGFILE select PROGFILE else "/dev/null";
procedure Say(s)
    fh := Open(PROG, "a");
    fprintf fh, "%o  %o\n", Cputime(), s;
    delete fh;
end procedure;

D := StringToInteger(DISC); N := StringToInteger(LEVEL);
PREC := StringToInteger(PR); TARGET := StringToInteger(TGT);
qmul := StringToInteger(QMUL); BOXSZ := StringToInteger(BOX);
NPROBE := StringToInteger(assigned NPR select NPR else "2500");
KMAX := 40;
DN := D*N; CC := ComplexField(PREC); RR := RealField(PREC); R4 := RealField(4);
W := { Integers() | StringToInteger(x) : x in Split(WSET, ",") };
Include(~W, 1);
Wstar := { Integers() | d : d in Divisors(DN) | GCD(d, DN div d) eq 1 };
error if not (W subset Wstar), "W is not a set of Hall divisors";
error if #Wstar ne 2*#W, Sprintf("W has index %o in the full group, not 2", #Wstar / #W);

t0 := Cputime();
A := QuaternionAlgebra(D); O := QuaternionOrder(A, N);
G := FuchsianGroup(O); H := UpperHalfPlane();
_, DD, FD := RecentredDisc(G, H, PREC, 3);
rhoG := Max([Abs(ComplexValue(v)) : v in FD]);
pt := ComplexValue(DiscToPlane(H, DD ! 0) : Precision := PREC);
mrep := FuchsianMatrixRepresentation(A : Precision := PREC);
tdom := Cputime(t0);

Xw := CreateShimuraQuot(D, N, W);     gw, esw := QuotientSignature(Xw);
Xs := CreateShimuraQuot(D, N, Wstar); gs, ess := QuotientSignature(Xs);
printf "=== X_0^%o(%o) / %o ===\n", D, N, Sort(SetToSequence(W));
printf "quotient signature (%o; %o) genus %o | star signature (%o; %o) genus %o\n",
    gw, esw, gw, gs, ess, gs;
printf "vol %o rho(Gamma) %o sides %o prec %o (domain %o s)\n",
    R4!ArithmeticVolume(G), R4!rhoG, #FD, PREC, R4!tdom;
Say(Sprintf("domain built, rho(Gamma) = %o", R4!rhoG));
error if gs ne 0, "the star curve does not have genus 0, so there is no hauptmodul to write over";

// ---- Atkin-Lehner elements and the normaliser's radius ----
t1 := Cputime();
Ws := SidePairingMatrices(G, DD, mrep, pt, CC);
Als := []; almat := AssociativeArray();
for m in Sort([ d : d in Wstar | d gt 1 ]) do
    ok, x := ALElement(O, m, BOXSZ);
    error if not ok, Sprintf("no Atkin-Lehner element of reduced norm %o", m);
    M := DiscMatrix(mrep, x, CC ! pt, CC);
    Append(~Als, <m, M>); almat[m] := M;
end for;
rmax := RR ! 0;
for i in [1..NPROBE] do
    r := (RR ! (0.05 + 0.95*((i*5) mod 23)/23)) * rhoG;
    w := CC ! (r * Exp(2*Pi(RR)*CC.1*i*(Sqrt(RR ! 5)-1)/2));
    _, v := ReduceOverStar(Ws, Als, w, CC);
    if Abs(v) gt rmax then rmax := Abs(v); end if;
end for;
rad := Min(RR ! 0.9995, 1.03*rmax);
NT := TermsFor(rad, TARGET);
printf "normaliser: rho(Gamma*) %o (probe %o points), sampling radius %o, NT %o (against %o on Gamma) (%o s)\n",
    R4!rmax, NPROBE, R4!rad, NT, TermsFor(rhoG, TARGET), R4!Cputime(t1);
Say(Sprintf("rho(Gamma*) = %o, NT = %o", R4!rmax, NT));

// ---- samples, labelled by the coset of the normaliser ----
t2 := Cputime();
samp := StarSamples(Ws, Als, rad, qmul*2*NT, CC);
for am in Als do
    samp cat:= StarSamplesVia(Ws, Als, am[2], am[1], rad, (qmul*NT) div 2, CC);
end for;
labs := AssociativeArray();
for s in samp do
    if IsDefined(labs, s[4]) then labs[s[4]] +:= 1; else labs[s[4]] := 1; end if;
end for;
printf "samples: %o, by coset %o (%o s)\n", #samp,
    [ <m, IsDefined(labs, m) select labs[m] else 0> : m in Sort(SetToSequence(Wstar)) ],
    R4!Cputime(t2);
Say(Sprintf("%o samples", #samp));
nW := &+([0] cat [ labs[m] : m in Wstar | m in W and IsDefined(labs, m) ]);
error if nW lt 2*NT, Sprintf("only %o samples in the subgroup's own coset, need %o", nW, 2*NT);

sgnW := AssociativeArray(); sgnS := AssociativeArray(); sgnM := AssociativeArray();
for m in Wstar do
    sgnS[m] := 1;
    if m in W then sgnW[m] := 1; sgnM[m] := 1; else sgnM[m] := -1; end if;
end for;
NR := 2*NT;

// ---- the hauptmodul of the star ----
function DiscPoint(i, f)
    return CC ! ((RR ! (0.12 + 0.45*((i*3) mod 7)/7) * f * rad) *
                 Exp(2*Pi(RR)*CC.1*i*(Sqrt(RR ! 5)-1)/2));
end function;

KT := 0; dT := 0;
for k in [2..KMAX by 2] do
    d := DimMkSig(k, gs, ess);
    if d ge 2 then KT := k; dT := d; break; end if;
end for;
error if KT eq 0, "no weight below the cap gives a star space of dimension at least two";
t3 := Cputime();
Mt := SystemLabelled(samp, KT, NT, sgnS, NR, CC);
if dT gt 2 then Mt := VerticalJoin(Mt, VanishRows([ DiscPoint(100+i, 1) : i in [1..dT-2] ], NT, CC)); end if;
Say(Sprintf("star pencil %o x %o; decomposing", Nrows(Mt), Ncols(Mt)));
bt, svt, gapt := KernelBasis(Mt, 2, CC);
Say(Sprintf("hauptmodul done, gap %o", R4!gapt));
printf "hauptmodul: weight %o, star dim %o, %o zeros imposed, rows %o, gap %o, kernel %o %o | next %o (%o s)\n",
    KT, dT, dT-2, Nrows(Mt), R4!gapt, R4!svt[1], R4!svt[2], R4!svt[3], R4!Cputime(t3);
error if gapt lt 1000, "the pencil on the star did not resolve; the hauptmodul would be meaningless";

// ---- y: the smallest weight carrying both an anti-invariant form and a star form ----
KY := 0; dminus := 0;
for k in [2..KMAX by 2] do
    dm := DimMkSig(k, gw, esw) - DimMkSig(k, gs, ess);
    if dm ge 1 and DimMkSig(k, gs, ess) ge 1 then KY := k; dminus := dm; break; end if;
end for;
error if KY eq 0, "no weight below the cap carries an anti-invariant form over a nonzero star form";
t4 := Cputime();
My := SystemLabelled(samp, KY, NT, sgnM, NR, CC);
Say(Sprintf("anti-invariant %o x %o; decomposing", Nrows(My), Ncols(My)));
bym, svy, gapy := KernelBasis(My, dminus, CC);
Say(Sprintf("anti-invariant done, gap %o", R4!gapy));
dplus := DimMkSig(KY, gs, ess);
if KY eq KT then
    bh := [ bt[1] ]; gaph := gapt;
else
    Mh := SystemLabelled(samp, KY, NT, sgnS, NR, CC);
    Say(Sprintf("star at the y weight %o x %o; decomposing", Nrows(Mh), Ncols(Mh)));
    bh, svh, gaph := KernelBasis(Mh, dplus, CC);
end if;
Say(Sprintf("y stage done, gaps %o %o", R4!gapy, R4!gaph));
printf "y: weight %o, anti-invariant dim %o (gap %o), invariant dim %o (gap %o) (%o s)\n",
    KY, dminus, R4!gapy, dplus, R4!gaph, R4!Cputime(t4);
error if gapy lt 1000 or gaph lt 1000, "the weight-KY spaces did not resolve";

// ---- values of t and y at points of the disc ----
// y^2 = P(t)/Q(t) holds identically on the disc, so the points need not lie in any fundamental
// domain; they are kept well inside the sampling circle, where the truncated series is far more
// accurate than at its edge.
npts := 8*(2*gw + 8);
ts := []; ys := [];
for i in [1..npts] do
    w := DiscPoint(i, 1);
    d2 := Ser(bt[2], w); hv := Ser(bh[1], w);
    if Abs(d2) lt 1e-30 or Abs(hv) lt 1e-30 then continue; end if;
    Append(~ts, Ser(bt[1], w)/d2);
    Append(~ys, Ser(bym[1], w)/hv);
end for;
printf "evaluation: %o points, |t| from %o to %o\n", #ts,
    R4!Min([Abs(x) : x in ts]), R4!Max([Abs(x) : x in ts]);

// ---- fit y^2 = P(t)/Q(t), raising the degrees until a one-dimensional kernel appears ----
P := []; Qd := []; dp := 0; dq := 0; res := 1.0; nxt := 0.0; found := false;
for dd in [1..2*gw+6] do
    for pq in [[dd,dd],[dd+1,dd],[dd,dd+1]] do
        if pq[1]+pq[2]+2 gt #ts then continue; end if;
        a, b, s1, s2 := FitRational(ts, ys, pq[1], pq[2], CC);
        if s1 lt 1e-12 and s2 gt 1e-4 then
            P := a; Qd := b; dp := pq[1]; dq := pq[2]; res := s1; nxt := s2; found := true; break;
        end if;
    end for;
    if found then break; end if;
end for;
Say(found select Sprintf("fit found, degrees %o %o", dp, dq) else "fit FAILED");
error if not found, "no rational function of the hauptmodul fitted y^2";
printf "fit: y^2 = P(t)/Q(t) with deg P <= %o, deg Q <= %o | residual %o, next singular value %o\n",
    dp, dq, R4!res, R4!nxt;

// ---- branch points: the roots of P*Q of ODD multiplicity, plus infinity when the degrees
// differ by an odd number.  A genus-g double cover of the line has exactly 2g+2 of them.
PC<T> := PolynomialRing(CC);
fp := &+[ P[i]*T^(i-1) : i in [1..#P] ];
fq := &+[ Qd[i]*T^(i-1) : i in [1..#Qd] ];
cls, lead := ClusteredRoots(Coefficients(fp*fq), 1e-8, CC);
br := [ c[1] : c in cls | IsOdd(c[2]) ];
atinf := IsOdd(Degree(fp) - Degree(fq));
printf "root clusters of P*Q: %o\n", [ <ComplexField(8)!c[1], c[2]> : c in cls ];
printf "BRANCHCOUNT %o expected %o%o\n", #br + (atinf select 1 else 0), 2*gw+2,
    atinf select " (one at infinity)" else "";
error if #br + (atinf select 1 else 0) ne 2*gw+2, "wrong number of branch points";

outdir := "/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/vvdata/weyl-campaign/voight-willis/out";
outf := Sprintf("%o/branch_%o_%o_%o.m", outdir, D, N, &cat[Sprint(m) cat "-" : m in Sort(SetToSequence(W))]);
fh := Open(outf, "w");
fprintf fh, "// branch points of X_0^%o(%o)/%o over the star, from the expansions\n",
    D, N, Sort(SetToSequence(W));
fprintf fh, "// weights %o and %o, NT %o, sampling radius %o, working precision %o\n",
    KT, KY, NT, rad, PREC;
fprintf fh, "D := %o; N := %o; W := %o; genus := %o; atinf := %o; prec := %o;\n",
    D, N, Sort(SetToSequence(W)), gw, atinf select "true" else "false", PREC;
fprintf fh, "CCb := ComplexField(%o);\nbranch := [CCb|\n", PREC;
for i in [1..#br] do
    fprintf fh, "  CCb ! [%o, %o]%o\n", Real(br[i]), Imaginary(br[i]), i lt #br select "," else "";
end for;
fprintf fh, "];\n";
delete fh;
printf "branch points written to %o\n", outf;
Say("branch points written");

// ---- comparison with an equation we already hold, or a published one, as SETS up to Mobius ----
function Finitise(zs, hasinf, c, CC)
    out := [ CC | 1/(z - c) : z in zs ];
    if hasinf then Append(~out, CC ! 0); end if;
    return out;
end function;

refs := [* *];
fn := Sprintf("/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/data/models/models_%o_%o.m", D, N);
if FileExists(fn) then
    models := eval (Read(fn) cat "\nreturn models;");
    key := Sort(SetToSequence(W));
    if IsDefined(models, key) and #models[key] ne 0 and Type(models[key][1][2]) ne MonStgElt then
        e := models[key][1];
        refs cat:= [* <"stored model", e[2] + (#e ge 3 select e[3]^2/4 else 0)> *];
    end if;
end if;
if assigned PUBF then
    PQ<x> := PolynomialRing(Rationals());
    refs cat:= [* <"published", eval PUBF> *];
end if;

cshift := CC ! (RR ! 0.3718 + 0.5171*CC.1);
mine := Finitise(br, atinf, cshift, CC);
for rf in refs do
    rts := [ CC | z[1] : z in Roots(ChangeRing(rf[2], CC)) ];
    hinf := IsOdd(Degree(rf[2]));
    if #rts + (hinf select 1 else 0) ne #mine then
        printf "MOBIUS %o: branch set has %o points against %o -- not comparable\n",
            rf[1], #rts + (hinf select 1 else 0), #mine;
        continue;
    end if;
    ok, resid := MobiusMatch(mine, Finitise(rts, hinf, cshift, CC), CC);
    printf "MOBIUS %o: %o (total residual %o over %o points, %o free after the three that fix the map)\n",
        rf[1], ok select "MATCH" else "NO MATCH", R4!resid, #mine, #mine-3;
end for;
printf "MODEL DONE D=%o N=%o W=%o\n", D, N, Sort(SetToSequence(W));
quit;
