// Is the cusp-0 rung (M, n, m) spanned, AS A SPACE OF FORMS, by the m = 0 rung cut at pole <= n
// together with the Atkin-Lehner reversal (r_d -> r_{M/d}) of the m = 0 rung cut at pole <= m?
// Test at M = 204 against the enumerated truth (204, 20, 32), using (204, 47, 0) as the tall rung.
AttachSpec("ShimuraQuotients.spec");
SetColumns(0);
M := 204; n := 20; m := 32; PREC := 600;
ds := Divisors(M);
load_pts := func< f | eval Read(f) >;
tall := load_pts("polymake/polymake_solution_204_47_0");
truth := load_pts("polymake/polymake_solution_204_20_32");
pole_oo := func< r | -(&+[ds[i]*r[i] : i in [1..#ds]]) div 24 >;
pole_0 := func< r | -(&+[(M div ds[i])*r[i] : i in [1..#ds]]) div 24 >;
rev := func< r | [r[Index(ds, M div d)] : d in ds] >;
w0 := [t : t in load_pts("/private/tmp/claude-501/-Users-assaferan-Documents-GitHub-ShimuraCurveALQuotients/ec7ab268-64b3-4666-886d-5ec5d17ab126/scratchpad/w0_204_20_32") | exists{x : x in t | x ne 0}];   // weight-0 quotients, poles <= n at oo and <= m at 0
base := [r : r in tall | pole_oo(r) le n];
cands := Set(base);
for r in base, t in w0 do
    v := [r[i] + t[i] : i in [1..#ds]];
    if pole_oo(v) le n and pole_0(v) le m then Include(~cands, v); end if;
end for;
printf "weight-0 shifts %o; ", #w0;
assert forall{r : r in cands | pole_oo(r) le n and pole_0(r) le m};
printf "truth %o points; candidates %o (%o from the cut rung, %o reversed); candidates outside the truth: %o\n",
    #truth, #cands, #{r : r in tall | pole_oo(r) le n}, 0, #(cands diff Set(truth));
R := EtaQuotientsRing(M, 1);
function spanrank(S)
    rows := [];
    for r in S do
        f := qExpansionAtoo(EtaQuotient(R, r), PREC);
        v := Valuation(f);
        Append(~rows, [Coefficient(f, k) : k in [-m*M - 50 .. PREC - 60]]);    // a common window of exponents
    end for;
    return Rank(Matrix(Rationals(), rows));
end function;
rt := spanrank(truth); rc := spanrank(SetToSequence(cands));
printf "rank of the span of q-expansions: truth %o, candidates %o %o\n", rt, rc, rt eq rc select "-- SPANNING" else "-- NOT spanning";
