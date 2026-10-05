// The fixed-m ladder: is the enumerated cusp-0 rung (M, nhi, m) spanned, as forms, by the lower
// enumerated rung (M, nlo, m) shifted (up to J times) by the weight-0 quotients with poles at oo only?
// magma -b M:=204 NLO:=20 NHI:=45 MM:=32 W0:=file ladder.m
AttachSpec("ShimuraQuotients.spec");
SetColumns(0);
M := StringToInteger(M); nlo := StringToInteger(NLO); nhi := StringToInteger(NHI); m := StringToInteger(MM);
ds := Divisors(M);
load_pts := func< f | eval Read(f) >;
low := load_pts(Sprintf("polymake/polymake_solution_%o_%o_%o", M, nlo, m));
truth := load_pts(Sprintf("polymake/polymake_solution_%o_%o_%o", M, nhi, m));
w0 := [t : t in load_pts(W0) | exists{x : x in t | x ne 0}];
pole_oo := func< r | -(&+[ds[i]*r[i] : i in [1..#ds]]) div 24 >;
pole_0 := func< r | -(&+[(M div ds[i])*r[i] : i in [1..#ds]]) div 24 >;
cands := Set(low); layer := Set(low);
for j in [1..3] do
    nxt := {};
    for r in layer, t in w0 do v := [r[i] + t[i] : i in [1..#ds]]; if pole_oo(v) le nhi then Include(~nxt, v); end if; end for;
    layer := nxt; cands join:= layer;
end for;
assert forall{r : r in cands | pole_oo(r) le nhi and pole_0(r) le m};
R := EtaQuotientsRing(M, 1);
PREC := 2*(nhi + m) + 400;
function spanrank(S)
    rows := [[Coefficient(f, k) : k in [-nhi .. PREC - nhi - 10]] where f := qExpansionAtoo(EtaQuotient(R, r), PREC) : r in S];
    return Rank(Matrix(Rationals(), rows));
end function;
rt := spanrank(truth); rc := spanrank(SetToSequence(cands));
printf "M=%o m=%o: rung %o (%o pts) -> rung %o (truth %o pts) via %o oo-shifts: %o candidates, %o outside the truth; span rank truth %o, candidates %o %o\n",
    M, m, nlo, #low, nhi, #truth, #w0, #cands, #(cands diff Set(truth)), rt, rc, rt eq rc select "SPANNING" else "NOT spanning";
