// Dual-graph filter (DualGraph.m).  Spec: docs/superpowers/specs/2026-09-29-dual-graph-filter-design.md.
// Why each number is trustworthy: genus identities are b_1(G_W) = GenusShimuraCurveQuotient (an
// independent formula); the K4 and the controls reproduce the 2026-09-29 pilot
// (handoff_2026-09-29/graph-pilot/); the Stankewicz data is raw output of fiber() from
// github.com/fsaia/GenusAtMost2.  Labels of LeftIdealClasses are not deterministic, so two builds
// are compared only up to isomorphism, never with eq.

// A deliberately wrong copy of a data tuple.  "noreverse" and "wrongreversal" keep the origin map
// equivariant, so their effect does not depend on the (non-deterministic) order of LeftIdealClasses;
// "trivialvertexAL" breaks equivariance, so its genus-failure count varies with the labelling
// (279-305 of 822 observed) and only "> 0" is asserted.
function dg_corrupt(data, p, variant)
    d := data;
    if variant eq "noreverse" then            // w_p trivial on edges: terminus = org[e]
        d[4] := [x[1] eq p select <p, [1..d[2]]> else x : x in d[4]];
    elif variant eq "wrongreversal" then      // w_p on edges replaced by w_q0, q0 the least other prime
        q0 := Min([x[1] : x in d[4] | x[1] ne p]);
        d[4] := [x[1] eq p select <p, [y[2] : y in d[4] | y[1] eq q0][1]> else x : x in d[4]];
    else                                      // "trivialvertexAL": a w_q that moves vertices made trivial on them
        q := [x[1] : x in d[5] | x[2] ne [1..d[1]]][1];
        d[5] := [x[1] eq q select <q, [1..d[1]]> else x : x in d[5]];
    end if;
    return d;
end function;

// fixed-point count of every element of the AL group, from a list of <q, images>: invariant under
// relabelling, so it compares builds (and Magma's Brandt modules) whose bases cannot be aligned
function dg_fixcounts(list, n)
    out := [];
    for S in Subsets({x[1] : x in list}) do
        img := [1..n];
        for q in S do qi := [x[2] : x in list | x[1] eq q][1]; img := [qi[i] : i in img]; end for;
        Append(~out, <&*[Integers() | q : q in S], #[i : i in [1..n] | img[i] eq i]>);
    end for;
    return Sort(out);
end function;

procedure test_DualGraphData()
    printf "Testing DualGraphData...";
    // (D, N, p, h, h'): class numbers at levels N and Np; h' - 2h + 1 = genus of X_0^D(N)
    for c in [<6,35,2,8,24>, <6,35,3,4,16>, <10,33,5,4,24>, <14,5,7,1,4>] do
        data := DualGraphData(c[1], c[2], c[3] : CacheDir := "none");
        assert data[1] eq c[4] and data[2] eq c[5];
    end for;
    // the full-scan path (the fallback of the class lookup) agrees with the keyed path
    for c in [<6,35,2>, <10,33,5>] do
        a := DualGraphData(c[1], c[2], c[3] : CacheDir := "none");
        b := DualGraphData(c[1], c[2], c[3] : ForceScan := true);
        assert dg_fixcounts(a[4], a[2]) eq dg_fixcounts(b[4], b[2]);
        assert dg_fixcounts(a[5], a[1]) eq dg_fixcounts(b[5], b[1]);
    end for;
    // validation rejects a vertex action that the origin map does not intertwine
    bad := dg_corrupt(DualGraphData(6, 35, 2 : CacheDir := "none"), 2, "trivialvertexAL");
    rejected := false;
    try ValidateDualGraphData(6, 35, 2, bad); catch e rejected := true; end try;
    assert rejected;
    // not applicable: D = 1, N not squarefree
    assert not DualGraphApplicable(1, 97) and not DualGraphApplicable(6, 25) and DualGraphApplicable(6, 35);
    printf "Done!\n";
end procedure;
test_DualGraphData();

procedure test_DualGraphDiskCache()
    printf "Testing DualGraph disk cache...";
    dir := Sprintf("/tmp/dualgraph_test_%o/nested", Getpid());
    System(Sprintf("rm -rf '/tmp/dualgraph_test_%o'", Getpid()));
    ClearDualGraphCache();
    d1 := DualGraphData(6, 35, 2 : CacheDir := dir);
    path := dir cat "/dualgraph_v1_6_35_2";
    assert OpenTest(path, "r");
    assert Pipe(Sprintf("ls '%o' | grep -c tmp || true", dir), "") eq "0\n";     // no temp file left
    ClearDualGraphCache();
    d2 := DualGraphData(6, 35, 2 : CacheDir := dir);
    assert d1 eq d2;                                   // read back from disk, same labelling
    // a corrupt file is rejected on read
    bad := dg_corrupt(d1, 2, "trivialvertexAL");
    Write(path, Sprint(bad, "Magma") : Overwrite := true);
    ClearDualGraphCache();
    rejected := false;
    try _ := DualGraphData(6, 35, 2 : CacheDir := dir); catch e rejected := true; end try;
    assert rejected;
    // a file with another format stamp is never read
    System(Sprintf("rm -f '%o'", path));
    Write(dir cat "/dualgraph_v0_6_35_2", "this is not Magma" : Overwrite := true);
    ClearDualGraphCache();
    _ := DualGraphData(6, 35, 2 : CacheDir := dir);
    assert OpenTest(path, "r");
    System(Sprintf("rm -rf '/tmp/dualgraph_test_%o'", Getpid()));
    ClearDualGraphCache();
    printf "Done!\n";
end procedure;
test_DualGraphDiskCache();

// simple graph from a multigraph, subdividing every edge (a loop by two points), for IsIsomorphic
function dg_subdivided(nv, edges)
    n := nv; E := {};
    for e in edges do
        if e[1] eq e[2] then E join:= {{e[1], n+1}, {n+1, n+2}, {n+2, e[1]}}; n +:= 2;
        else E join:= {{e[1], n+1}, {n+1, e[2]}}; n +:= 1; end if;
    end for;
    return Graph< n | E >;
end function;

// two dual-graph data tuples at p give isomorphic quotient graphs for every AL subgroup
function dg_same_graphs(a, b, D, N, p)
    for W in ALSubgroups(D*N) do
        n1, e1 := DualGraphQuotientFromData(a, p, W[1]);
        n2, e2 := DualGraphQuotientFromData(b, p, W[1]);
        if not IsIsomorphic(dg_subdivided(n1, e1), dg_subdivided(n2, e2)) then return false; end if;
    end for;
    return true;
end function;

// (D, N, W, g) rows; read at top level: eval inside a function or procedure of an eval'd test crashes Magma 2.29-4
dg_genus_rows := eval Read("tests/dualgraph_genus_data.txt");

procedure test_DualGraphQuotient(rows)
    printf "Testing DualGraphQuotientFromData...";
    // genus identities b_1(G_W) = g, D*N <= 700 (822 of them; all 988 in tests/_offline/DualGraphFull.m)
    n := 0;
    for r in [r : r in rows | r[1]*r[2] le 700] do
        for p in PrimeDivisors(r[1]) do
            nv, edges := DualGraphQuotientFromData(DualGraphData(r[1], r[2], p : CacheDir := "none"), p, r[3]);
            assert #edges - nv + 1 eq r[4] or (#edges eq 0 and r[4] eq 0);
            n +:= 1;
        end for;
    end for;
    assert n eq 822;
    // the full-scan build gives the same graphs, for every W
    for c in [<6,35,2>, <10,33,5>] do
        assert dg_same_graphs(DualGraphData(c[1], c[2], c[3] : CacheDir := "none"),
                              DualGraphData(c[1], c[2], c[3] : ForceScan := true), c[1], c[2], c[3]);
    end for;
    // negative controls: each deliberately wrong graph (dg_corrupt) breaks at least one identity
    for variant in ["noreverse", "wrongreversal", "trivialvertexAL"] do
        fails := 0;
        for r in [r : r in rows | r[1]*r[2] le 700] do
            for p in PrimeDivisors(r[1]) do
                data := DualGraphData(r[1], r[2], p : CacheDir := "none");
                if variant eq "trivialvertexAL" and forall{x : x in data[5] | x[2] eq [1..data[1]]} then continue; end if;
                nv, edges := DualGraphQuotientFromData(dg_corrupt(data, p, variant), p, r[3]);
                b1 := #edges eq 0 select 0 else #edges - nv + 1;
                if b1 ne r[4] then fails +:= 1; end if;
            end for;
        end for;
        printf " %o:%o", variant, fails;
        assert fails gt 0;
    end for;
    printf " Done!\n";
end procedure;
test_DualGraphQuotient(dg_genus_rows);
