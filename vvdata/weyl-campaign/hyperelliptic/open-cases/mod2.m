// Mod-2 Galois image test.  For a hyperelliptic curve of genus g over Q, J[2] is the even subsets
// of the 2g+2 Weierstrass points modulo the complement, so Frobenius acts through a permutation
// and its characteristic polynomial on J[2], the Weil polynomial mod 2, is prod (x^{c_i}+1)/(x+1)^2
// for a partition (c_i) of 2g+2.  A Weil polynomial mod 2 outside that set proves X not hyperelliptic.
AttachSpec("ShimuraQuotients.spec"); SetColumns(0);
R<x> := PolynomialRing(GF(2));
allowed := AssociativeArray();
for g in [3, 4, 5, 6] do allowed[g] := {(&*[x^c + 1 : c in P]) div (x+1)^2 : P in Partitions(2*g + 2)}; end for;
cs := eval Read("data/curves_after_UpdateCurves8.dat");
open := [X : X in cs | not assigned X`IsSubhyp];
nr := 0;
for X in open do
    g := X`g; bound := g le 3 select 60 else 30;
    bad := []; seen := [];
    for p in PrimesUpTo(bound) do
        if (X`D * X`N) mod p eq 0 then continue; end if;
        f := R ! WeilPolynomial(X, p);
        Append(~seen, <p, f>);
        if not (f in allowed[g]) then Append(~bad, <p, f>); end if;
    end for;
    if #bad gt 0 then nr +:= 1; end if;
    printf "X_0(%o,%o)/%o g=%o: %o primes, %o forbidden%o\n", X`D, X`N, Sort(Setseq(X`W)), g, #seen, #bad, #bad gt 0 select Sprintf(" => NOT hyperelliptic, e.g. p=%o: %o", bad[1][1], bad[1][2]) else "";
end for;
printf "MOD2 DONE: %o of %o ruled out\n", nr, #open;
quit;
