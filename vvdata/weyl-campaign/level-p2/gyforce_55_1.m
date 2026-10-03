AttachSpec("ShimuraQuotients.spec");
SetColumns(0);
SetVerbose("ShimuraQuotients", 1);
D := 55; N := 1; extra := [-27];
gy := [
<-3, -3/2>,
<-12, Infinity()>,
<-15, 1>,
<-27, 0>,
<-60, -1>,
<-67, 8/3>,
<-88, -2>,
<-115, -4>,
<-163, -5/6>,
<-187, 1/2>,
<-235, -4/11>,
<-715, -4/3>
];
gyv := AssociativeArray(); for t in gy do gyv[t[1]] := t[2]; end for;
Xstar := CreateShimuraQuot(D, N, Set(Divisors(D*N)));
Xstar`g := GenusShimuraCurveQuotient(D, N, Xstar`W); Xstar`CurveID := 0;
curves := GetQuotientsAndGenera([Xstar]);
_ := exists(star){c : c in curves | IsStarCurve(c)};
fs := BorcherdsForms(star, curves : Prec := 100);
d_divs := &cat[[T[1]: T in DivisorOfBorcherdsForm(f, star)] : f in [fs[-1], fs[-2]]];
must_use := Set(d_divs) join {t[1] : t in gy};
cm := CandidateDiscriminants(star, curves : Keep := must_use);
rat := cm[1]; quad := cm[2];
for d in extra do
    d0 := FundamentalDiscriminant(d); _, f := IsSquare(d div d0);
    OK := MaximalOrder(QuadraticField(d)); O := sub<OK | f>;
    npts := NumberOfOptimalEmbeddings(O, D, N) div 2^#PrimeDivisors(D*N);
    printf "forcing d = %o (conductor %o): %o star point(s)\n", d, f, npts;
    if not exists{u : u in rat cat quad | u[1] eq d} then
        if npts eq 1 then Append(~rat, <d, f, 1>); else Append(~quad, <d, f, 2>); end if;
    end if;
end for;
abs_tab, all_cm_pts := AbsoluteValuesAtCMPoints(star, curves, [rat, quad], fs : MaxNum := 12, Prec := 100, Exclude := {}, Include := must_use);
ReduceTable(abs_tab);
tab := ValuesAtCMPoints(abs_tab, all_cm_pts);
printf "KEYS %o sIndex %o sTildeIndex %o\n", tab`Keys_fs, tab`sIndex, tab`sTildeIndex;
discs := tab`Discs; srow := tab`Values[tab`sIndex]; strow := tab`Values[tab`sTildeIndex];
function mobius(z0, z1, z2, z)
    num := (z eq Infinity() or z0 eq Infinity()) select 1 else z - z0;
    den := (z eq Infinity() or z1 eq Infinity()) select 1 else z - z1;
    c1  := (z2 eq Infinity() or z1 eq Infinity()) select 1 else z2 - z1;
    c2  := (z2 eq Infinity() or z0 eq Infinity()) select 1 else z2 - z0;
    if den*c2 eq 0 then return Infinity(); end if;
    return (num*c1)/(den*c2);
end function;
idx := AssociativeArray(); for i->d in discs do idx[d] := i; end for;
have := [t : t in gy | IsDefined(idx, t[1])];
ref := []; for t in have do if #ref eq 3 then break; end if; if forall{r : r in ref | r[2] ne t[2]} and t[1] notin extra then Append(~ref, t); end if; end for;
for name in ["s", "stilde"] do
    row := name eq "s" select srow else strow;
    z := [row[idx[t[1]]] : t in ref]; w := [t[2] : t in ref];
    for t in have do
        if t[1] in {r[1] : r in ref} then continue; end if;
        got := mobius(z[1], z[2], z[3], row[idx[t[1]]]); want := mobius(w[1], w[2], w[3], t[2]);
        printf "CHECK %o d %-6o %o\n", name, t[1], got eq want select "OK" else Sprintf("MISMATCH got %o want %o", got, want);
    end for;
end for;
