// Does sol(74) + t-shifts span the same space as the full sol(128)?
// Expected dim of {pole <= 128 at oo} = 108 (rank(sol_74) = 54, +k per rung).
AttachSpec("ShimuraQuotients.spec");
Ldata := ShimuraCurveLattice(77, 2);
disc := #(Ldata`Ldual/Ldata`L);
M := 308; Prec := 700;
R := EtaQuotientsRing(M, disc);
ds := Divisors(M);
sol74 := Sort(eval Read("polymake/polymake_solution_308_74_0"));
W0 := [ r : r in eval Read(W0F) | r ne [ 0 : d in ds ] ];
printf "sol74: %o, weight-0 shifts: %o\n", #sol74, #W0;
poles := [ -(&+[ds[i]*t[i] : i in [1..#ds]])/24 : t in W0 ];
printf "shift poles: %o\n", poles;
cand := {@ r : r in sol74 @};
for t in W0 do
    for r in sol74 do Include(~cand, [ r[i] + t[i] : i in [1..#ds] ]); end for;
end for;
printf "candidate pool: %o\n", #cand;
etas := [ EtaQuotient(R, r) : r in cand ];
qexps := [ qExpansionAtoo(eta, Prec) : eta in etas ];
_<q> := Universe(qexps);
min_v := Minimum([Valuation(f) : f in qexps]);
coeffs := Matrix(Rationals(), [AbsEltseq(q^(-min_v)*f : FixedLength) : f in qexps]);
rk := Rank(coeffs);
printf "rank of shifted union = %o (target dim 108)\n", rk;
quit;
