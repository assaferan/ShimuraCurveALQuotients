// Rank-preserving core of sol(420,145) via one transpose-echelon:
// pivot columns of A^T index a maximal independent subset of A's rows.
AttachSpec("ShimuraQuotients.spec");
Ldata := ShimuraCurveLattice(210, 1);
disc := #(Ldata`Ldual/Ldata`L);
M := 420; Prec := 400;
R := EtaQuotientsRing(M, disc);
ds := Divisors(M);
sol := Sort(eval Read("polymake/polymake_solution_420_145_0"));
poleof := func< r | -(&+[ds[i]*r[i] : i in [1..#ds]]) div 24 >;
Sort(~sol, func< a, b | poleof(b) - poleof(a) >);
qexps := [ qExpansionAtoo(EtaQuotient(R, r), Prec) : r in sol ];
_<q> := Universe(qexps);
min_v := Minimum([Valuation(f) : f in qexps]);
A := Matrix(Rationals(), [AbsEltseq(q^(-min_v)*f : FixedLength) : f in qexps]);
printf "matrix %o x %o built\n", Nrows(A), Ncols(A);
Et := EchelonForm(Transpose(A));
rk := Rank(Et);
idx := [ PivotColumn(Et, i) : i in [1..rk] ];
printf "rank %o, core indices selected\n", rk;
core := [ sol[i] : i in idx ];
out := "[ PowerSequence(IntegerRing()) |\n";
for i->r in core do
    out cat:= "[ " cat Join([Sprint(x) : x in r], ", ") cat " ]" cat (i lt #core select ",\n" else "\n");
end for;
out cat:= "]\n";
Write("polymake/tshift_core_420.txt", out : Overwrite);
printf "core %o vectors -> polymake/tshift_core_420.txt\n", #core;
quit;
