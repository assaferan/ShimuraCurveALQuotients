// A prediction of the level-p^2 analysis, testable from the lattice alone: the number of isotropic
// cosets of L^v/L.  At squarefree level N the count is 2N-1 (tests/VectorValuedForm.m).  The
// p^2-scaled plane has 3p^2 - 2p isotropic cosets including 0, and the D-part is anisotropic, so
// for p^2 || N the prediction is 3p^2 - 2p, NOT 2N - 1.
AttachSpec("ShimuraQuotients.spec");
for b in [[6,5],[6,25],[6,49],[6,7],[15,2],[15,4]] do
    D := b[1]; N := b[2];
    Ld := ShimuraCurveLattice(D, N);
    Q := ChangeRing(Ld`Q, Rationals()); dn := Ld`denom;
    n := 0;
    for g in Ld`disc_grp do
        v := ChangeRing(g @@ Ld`to_disc, Rationals());
        r := (v*Q, v)/(2*dn^2);
        if r eq Floor(r) then n +:= 1; end if;
    end for;
    p := [q : q in PrimeDivisors(N)][1];
    e := Valuation(N, p);
    pred := (e eq 1) select 2*N - 1 else (e eq 2 select 3*p^2 - 2*p else 0);
    printf "%o_%o: |disc group| = %-7o isotropic cosets = %-5o prediction %-5o %o\n",
           D, N, #Ld`disc_grp, n, pred, (pred eq 0) select "(no prediction)" else (n eq pred select "MATCHES" else "MISMATCH");
end for;
