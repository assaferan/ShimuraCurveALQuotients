// THE DECISIVE TEST.  Guo-Yang Example 36: on X = X0^14(5), with s the Hauptmodul of X/W_{14,5}
// normalised by s(tau_-4) = oo, s(tau_-11) = 1, s(tau_-35) = 0, Schofer's formula gives
//        s(tau_-280) = 5/16.
// Why this discriminates.  N = 5 SPLITS in Q(sqrt -4) and Q(sqrt -11) -- both FIRING -- while
// 5 | 35 and 5 | 280, so tau_-35 and tau_-280 are NON-firing.  A correction of delta*log 5 applied
// only at firing discriminants therefore CHANGES the cross-ratio and moves s(tau_-280).  Unlike a
// model test, no row rescaling can absorb it.
//
// main applies NO m = 0 term at odd N (the IsEven(N) guard), so this run is the UNCORRECTED pipeline.
//   * if it returns 5/16, the uncorrected Schofer formula is already right for CM VALUES, and the
//     correction the model path needs is compensating something else in the pipeline;
//   * if it does not, the correction is real mathematics and we can read off what it must be.
AttachSpec("ShimuraQuotients.spec");

function mobius(z0, z1, z2, z)
    num := (z eq Infinity() or z0 eq Infinity()) select 1 else z - z0;
    den := (z eq Infinity() or z1 eq Infinity()) select 1 else z - z1;
    c1  := (z2 eq Infinity() or z1 eq Infinity()) select 1 else z2 - z1;
    c2  := (z2 eq Infinity() or z0 eq Infinity()) select 1 else z2 - z0;
    if den*c2 eq 0 then return Infinity(); end if;
    return (num*c1)/(den*c2);
end function;

D := 14; N := 5;
d0 := -35; d1 := -4; d2 := -11;          // -> 0, oo, 1  (GY: s(-35)=0, s(-4)=oo, s(-11)=1)
target := -280; expected := Rationals()!(5/16);
keep := {d0, d1, d2, target};
for d in keep do
    df := FundamentalDiscriminant(d);
    printf "d = %-6o d_fund = %-6o : %o\n", d, df,
           (df mod N eq 0) select "NON-firing (5 | d)" else
           Sprintf("FIRING, 5 is %o", KroneckerSymbol(df,N) eq 1 select "split" else "inert");
end for;

Xstar := CreateShimuraQuot(D, N, Set(Divisors(D*N)));
Xstar`g := GenusShimuraCurveQuotient(D, N, Xstar`W); Xstar`CurveID := 0;
curves := GetQuotientsAndGenera([Xstar]);
assert exists(star){c : c in curves | IsStarCurve(c)};
t0 := Cputime();
tab := ValuesAtCMPoints(star, curves : Keep := keep);
printf "ValuesAtCMPoints: %o s\n", Cputime(t0);
discs := tab`Discs;
srow := tab`Values[tab`sIndex];
idx := AssociativeArray(); for i->d in discs do idx[d] := i; end for;
printf "table discs: %o\n", discs;
for d in keep do
    if not IsDefined(idx, d) then printf "MISSING disc %o\n", d; end if;
end for;
if &and[IsDefined(idx, d) : d in keep] then
    got := mobius(srow[idx[d0]], srow[idx[d1]], srow[idx[d2]], srow[idx[target]]);
    printf "\ns(tau_%o) computed = %o    published = %o    %o\n",
           target, got, expected, (got eq expected) select "MATCH" else "*** MISMATCH ***";
    printf "raw row values: %o\n", [<d, srow[idx[d]]> : d in Sort(Setseq(keep))];
end if;
quit;
