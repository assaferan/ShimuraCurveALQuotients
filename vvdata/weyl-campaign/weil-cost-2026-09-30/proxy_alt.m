// Candidate Weil proxy: drop the Qmax factor.  Measured per-prime cost went like n^0.8 with
// n = p^g and was nearly independent of Qmax, W and the level (13029 p=31, n=9.2e5: 206 s;
// 1071 p=13, n=4.8e6: 560 s; 1071 p=17, n=2.4e7: 2132 s).  So cost ~ sum_p p^(0.8 g).
AttachSpec("ShimuraQuotients.spec");
curves := eval Read("data/curves_after_UpdateCurves7.dat");
ceil := AssociativeArray(); ceil[3]:=53; ceil[4]:=53; ceil[5]:=37; ceil[6]:=29; ceil[7]:=23; ceil[8]:=17;
function alt(X)
    g := X`g; DN := X`D*X`N; Qmax := Max(X`W);
    b := WeilClassNumberPrimeBound(Qmax, g); bb := WeilBudgetPrimeBound(Qmax, g);
    if bb lt b then b := bb; end if;
    if IsDefined(ceil, g) and ceil[g] lt b then b := ceil[g]; end if;
    return &+[RealField(6) | (RealField(6)!p)^(0.8*g) : p in PrimesUpTo(b) | DN mod p ne 0];
end function;
open := [c : c in curves | not assigned c`IsSubhyp and c`g ge 3];
pr := [<alt(c), c`CurveID, c`D, c`N, c`g, #c`W, Max(c`W)> : c in open];
Sort(~pr, func<a, b | b[1] - a[1]>);
ids := [t[2] : t in pr];
printf "alt proxy: 13029 -> %o (rank %o), 1071 -> %o (rank %o) of %o\n",
    alt(curves[13029]), Position(ids, 13029), alt(curves[1071]), Position(ids, 1071), #open;
printf "top 12 by alt proxy (proxy, id, D, N, g, #W, Qmax):\n";
for t in pr[1..12] do printf "  %o\n", t; end for;
printf "genus histogram of top 50: %o\n", {* t[5] : t in pr[1..50] *};
exit;
