// Hyperellipticity of every stored CRV curve by the canonical quadrics: among the g(g+1)/2 products
// of a basis of holomorphic differentials there are (g-1)(g-2)/2 linear relations when the curve is
// hyperelliptic (canonical image a rational normal curve) and (g-2)(g-3)/2 otherwise (Noether).
AttachSpec("ShimuraQuotients.spec");
SetColumns(0);
P<x> := PolynomialRing(Rationals());
files := Split(Pipe("ls data/models/models_*.m", ""), "\n");
nh := 0; nn := 0; nerr := 0;
for f in files do
    if #f eq 0 then continue; end if;
    b := Split(Split(f, "/")[3], ".")[1]; b := b[8..#b];
    models := eval (Read(f) cat "\nreturn models;");
    for k in Sort([k : k in Keys(models)]) do
        for e in models[k] do
            if Type(e[2]) ne MonStgElt then continue; end if;
            t0 := Cputime();
            try
                fs := CRVDoubleCovers(e[3]);
                K := FibreProductFunctionField(fs);
                g := Genus(K);
                if g le 2 then
                    printf "%o %o genus %o: hyperelliptic (genus <= 2)  %.1os\n", b, k, g, Cputime(t0); nh +:= 1; continue;
                end if;
                V, mp := SpaceOfHolomorphicDifferentials(K);
                ws := [mp(V.i) : i in [1..g]];
                phis := [ws[i] / ws[1] : i in [1..g]];
                prods := [phis[i]*phis[j] : j in [i..g], i in [1..g]];
                nrel := Dimension(Relations(prods, Rationals()));
                hyp := nrel eq (g-1)*(g-2) div 2;
                nonhyp := nrel eq (g-2)*(g-3) div 2;
                assert hyp or nonhyp;
                printf "%o %o genus %o: %o (%o quadrics)  %.1os\n", b, k, g, hyp select "HYPERELLIPTIC" else "not hyperelliptic", nrel, Cputime(t0);
                if hyp then nh +:= 1; else nn +:= 1; end if;
            catch err
                printf "%o %o: ERROR %o\n", b, k, err`Object; nerr +:= 1;
            end try;
        end for;
    end for;
end for;
printf "SWEEP DONE: %o hyperelliptic, %o not, %o errors\n", nh, nn, nerr;
quit;
