AttachSpec("ShimuraQuotients.spec"); SetColumns(0);
for c in [<39, 11, 8, "429b1", [-132, -52]>, <55, 7, 8, "385a1", [-1540, -20]>] do
    D, N, degf, lab, ds := Explode(c); E := EllipticCurve(lab);
    printf "== %o: signs %o\n", lab, [<v, QuaternionicAtkinLehnerSign(E, v, D, N)> : v in PrimeDivisors(D*N)];
    for d in ds do
        comps, h := HeegnerHeightComponents(E, D, N, degf, d);
        printf "d = %o (h = %o): %o\n", d, h, [<x[1], x[2], x[3], x[4], x[5], Sprintf("%.5o", x[6])> : x in comps];
    end for;
end for;
