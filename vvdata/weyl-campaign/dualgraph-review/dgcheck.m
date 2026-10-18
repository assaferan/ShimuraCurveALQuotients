// Corroboration of paper/dualgraph-nonsquarefree.tex at small non-squarefree levels.
// Run from the worktree root: magma -b dgcheck.m < /dev/null > dgcheck.out
AttachSpec("ShimuraQuotients.spec");
SetColumns(0);

// (1) The Atkin-Lehner ideal at q = l^k, k >= 2, built by dg_prime, is not l^(k/2) O and is
//     the reflection: check through the public interface by building the data and reading
//     the vertex action of w_q: it must differ from the identity whenever w_q moves classes.
//     Also record h, h' and the unit orders (edge stabilisers = half the unit group).
for t in [<6, 25>, <22, 9>, <15, 16>, <14, 9>] do
    D := t[1]; N := t[2]; DN := D*N;
    full := {Integers()| d : d in Divisors(DN) | GCD(d, DN div d) eq 1};
    subs := [W : W in Subsets(full) | 1 in W and forall{<a, b> : a, b in W | AtkinLehnerMul(a, b, DN) in W}];
    for p in PrimeDivisors(D) do
        tt := Cputime();
        data := DualGraphData(D, N, p : CacheDir := "none");   // ValidateDualGraphData runs inside
        printf "D=%o N=%o p=%o: h=%o h'=%o build+validate %.1os; unit orders V=%o E=%o\n",
            D, N, p, data[1], data[2], Cputime(tt), Seqset(data[6]), Seqset(data[7]);
        for x in data[5] do
            printf "   vertex action of w_%o: %o fixed classes of %o\n", x[1], #[i : i in [1..data[1]] | x[2][i] eq i], data[1];
        end for;
        for x in data[4] do
            printf "   edge action of w_%o: %o fixed classes of %o\n", x[1], #[i : i in [1..data[2]] | x[2][i] eq i], data[2];
        end for;
        bad := 0;
        for W in subs do
            nv, edges := DualGraphQuotientFromData(data, p, W);
            b1 := #edges eq 0 select 0 else #edges - nv + 1;
            g := GenusShimuraCurveQuotient(D, N, W);
            if b1 ne g then bad +:= 1; printf "   MISMATCH W=%o: b1=%o genus=%o\n", W, b1, g; end if;
        end for;
        printf "   genus identities: %o subgroups, %o mismatches\n", #subs, bad;
    end for;
end for;

// (2) Negative control for the prime-power Atkin-Lehner ideal: replace the vertex action of w_25
//     at (6, 25, p = 2) by the identity (what l^(k/2) O would induce) and show the genus
//     identities then fail for some W containing 25.
data := DualGraphData(6, 25, 2 : CacheDir := "none");
d := data;
d[5] := [x[1] eq 25 select <25, [1..d[1]]> else x : x in d[5]];
d[4] := [x[1] eq 25 select <25, [1..d[2]]> else x : x in d[4]];
DN := 150;
full := {Integers()| d0 : d0 in Divisors(DN) | GCD(d0, DN div d0) eq 1};
subs := [W : W in Subsets(full) | 1 in W and forall{<a, b> : a, b in W | AtkinLehnerMul(a, b, DN) in W}];
bad := 0;
for W in subs do
    nv, edges := DualGraphQuotientFromData(d, 2, W);
    b1 := #edges eq 0 select 0 else #edges - nv + 1;
    if b1 ne GenusShimuraCurveQuotient(6, 25, W) then bad +:= 1; end if;
end for;
printf "control (w_25 made trivial at 6_25, p=2): %o of %o genus identities fail\n", bad, #subs;
printf "DONE\n";
quit;
