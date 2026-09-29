// All 988 genus identities b_1(G_W) = g of the 2026-09-29 pilot list (14 levels, up to (1974, 1)).
// Offline: about 15 s, dominated by the fibers of D = 1974 and D = 546.
dg_full_rows := eval Read("tests/dualgraph_genus_data.txt");
n := 0;
for r in dg_full_rows do
    for p in PrimeDivisors(r[1]) do
        nv, edges := DualGraphQuotientFromData(DualGraphData(r[1], r[2], p : CacheDir := "none"), p, r[3]);
        assert #edges - nv + 1 eq r[4] or (#edges eq 0 and r[4] eq 0);
        n +:= 1;
    end for;
end for;
assert n eq 988;
printf "DualGraphFull: %o genus identities\n", n;
