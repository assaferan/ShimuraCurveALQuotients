// End-to-end: the pipeline's Hauptmodul rows s, s~ with the conductor-4 points -240, -48 in the
// table, compared with Guo-Yang Table 45 by cross-ratios (as tests/_offline/GuoYangCheck.m does).
// Variant A: current code.  Variant B (env FIBRE_* set): the m = 0-type terms at -240 and -48
// replaced by the fibre sum of fibresum.m.
AttachSpec("ShimuraQuotients.spec");
SetColumns(0);
D := 15; N := 2;
Xstar := CreateShimuraQuot(D, N, Set(Divisors(D*N)));
Xstar`g := GenusShimuraCurveQuotient(D, N, Xstar`W); Xstar`CurveID := 0;
curves := GetQuotientsAndGenera([Xstar]);
_ := exists(star){c : c in curves | IsStarCurve(c)};
fs := BorcherdsForms(star, curves : Prec := 100);
gy := [ <-7, 1/4>, <-12, Infinity()>, <-15, 5/4>, <-28, 9/4>, <-40, 0>, <-48, -1/4>, <-52, 1>, <-60, -1/12>,
<-88, 4>, <-120, 2>, <-132, -1>, <-148, 1/25>, <-168, 2/3>, <-228, -1/9>, <-232, 144/121>, <-240, -25/12>,
<-280, 10>, <-312, 2/25>, <-340, 9/17>, <-372, -31/9>, <-408, 68/25>, <-420, 5/3>, <-520, -8/121>,
<-660, -5/11>, <-708, -841/121>, <-760, 450/529>, <-840, 40/27> ];
gyv := AssociativeArray(); for t in gy do gyv[t[1]] := t[2]; end for;
d_divs := &cat[[T[1]: T in DivisorOfBorcherdsForm(f, star)] : f in [fs[-1], fs[-2]]];
must_use := Set(d_divs) join {-7, -15, -52, -28, -60, -240, -48};
cm := CandidateDiscriminants(star, curves : Keep := must_use);
rat := cm[1]; quad := cm[2];
for t in [<-28, 2, 1>, <-240, 4, 1>, <-48, 4, 1>] do
    if not exists{u : u in rat | u[1] eq t[1]} then Append(~rat, t); end if;
end for;
abs_tab, all_cm_pts := AbsoluteValuesAtCMPoints(star, curves, [rat, quad], fs : MaxNum := 12, Prec := 100, Exclude := {}, Include := must_use);
printf "KEYS %o  sIndex %o sTildeIndex %o\n", abs_tab`Keys_fs, abs_tab`sIndex, abs_tab`sTildeIndex;
printf "DISCS %o\n", abs_tab`Discs;
for i in [1..#abs_tab`Values] do printf "RAW row %o: %o\n", i, abs_tab`Values[i]; end for;
ReduceTable(abs_tab);
for i in [1..#abs_tab`Values] do printf "RED row %o: %o\n", i, abs_tab`Values[i]; end for;
try
    tab := ValuesAtCMPoints(abs_tab, all_cm_pts);
catch e
    printf "ValuesAtCMPoints failed: %o\n", e`Object; quit;
end try;
discs := tab`Discs;
srow := tab`Values[tab`sIndex]; strow := tab`Values[tab`sTildeIndex];
function mobius(z0, z1, z2, z)
    num := (z eq Infinity() or z0 eq Infinity()) select 1 else z - z0;
    den := (z eq Infinity() or z1 eq Infinity()) select 1 else z - z1;
    c1  := (z2 eq Infinity() or z1 eq Infinity()) select 1 else z2 - z1;
    c2  := (z2 eq Infinity() or z0 eq Infinity()) select 1 else z2 - z0;
    if den*c2 eq 0 then return Infinity(); end if;
    return (num*c1)/(den*c2);
end function;
// frame on three fundamental points
ref := [-7, -15, -52];
idx := AssociativeArray(); for i->d in discs do idx[d] := i; end for;
for name in ["s", "stilde"] do
    row := name eq "s" select srow else strow;
    z := [row[idx[d]] : d in ref]; w := [gyv[d] : d in ref];
    printf "ROW %o: frame %o -> %o\n", name, z, w;
    for d in discs do
        if not IsDefined(gyv, d) or d in ref then continue; end if;
        got := mobius(z[1], z[2], z[3], row[idx[d]]); want := mobius(w[1], w[2], w[3], gyv[d]);
        printf "CHECK %o d %-5o value %-12o got %-10o want %-10o %o\n", name, d, row[idx[d]], got, want, got eq want select "OK" else "MISMATCH";
    end for;
end for;
