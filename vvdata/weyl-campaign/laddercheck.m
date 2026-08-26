// Rank of the augmented-core t-shift ladder at J=1 and J=2 (M=420).
// Targets: dim(249) = 209, dim(353) = 313 (dim = n - 40 from sol145 rank 105).
AttachSpec("ShimuraQuotients.spec");
Ldata := ShimuraCurveLattice(210, 1);
disc := #(Ldata`Ldual/Ldata`L);
M := 420; Prec := 800;
R := EtaQuotientsRing(M, disc);
ds := Divisors(M);
core := eval Read("polymake/tshift_core_420.txt");
w0 := [ r : r in eval Read("polymake/tshift_w0_420.txt") | r ne [ 0 : d in ds ] ];
t := w0[1];
printf "core %o vectors, shift pole %o\n", #core, -(&+[ds[i]*t[i] : i in [1..#ds]]) div 24;
for J in [1, 2] do
    cand := {@ @};
    for j in [0..J] do
        for r in core do Include(~cand, [ r[i] + j*t[i] : i in [1..#ds] ]); end for;
    end for;
    qexps := [ qExpansionAtoo(EtaQuotient(R, r), Prec) : r in cand ];
    _<q> := Universe(qexps);
    min_v := Minimum([Valuation(f) : f in qexps]);
    coeffs := Matrix(Rationals(), [AbsEltseq(q^(-min_v)*f : FixedLength) : f in qexps]);
    printf "J=%o: %o candidates, rank = %o (target %o)\n", J, #cand, Rank(coeffs), 105 + 104*J;
end for;
quit;
