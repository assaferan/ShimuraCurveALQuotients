// Dual-graph filter (DualGraph.m).  Spec: docs/superpowers/specs/2026-09-29-dual-graph-filter-design.md.
// Why each number is trustworthy: genus identities are b_1(G_W) = GenusShimuraCurveQuotient (an
// independent formula); the K4 and the controls reproduce the 2026-09-29 pilot
// (handoff_2026-09-29/graph-pilot/); the Stankewicz data is raw output of fiber() from
// github.com/fsaia/GenusAtMost2.  Labels of LeftIdealClasses are not deterministic, so two builds
// are compared only up to isomorphism, never with eq.

// A deliberately wrong copy of a data tuple.  "noreverse", "wrongreversal" and "twistedreversal" keep
// the origin map equivariant, so their effect does not depend on the (non-deterministic) order of
// LeftIdealClasses;
// "trivialvertexAL" breaks equivariance, so its genus-failure count varies with the labelling
// (279-305 of 822 observed) and only "> 0" is asserted.
function dg_corrupt(data, p, variant)
    d := data;
    if variant eq "noreverse" then            // w_p trivial on edges: terminus = org[e]
        d[4] := [x[1] eq p select <p, [1..d[2]]> else x : x in d[4]];
    elif variant eq "wrongreversal" then      // w_p on edges replaced by w_q0, q0 the least other prime
        q0 := Min([x[1] : x in d[4] | x[1] ne p]);
        d[4] := [x[1] eq p select <p, [y[2] : y in d[4] | y[1] eq q0][1]> else x : x in d[4]];
    elif variant eq "twistedreversal" then    // w_p on edges replaced by w_p w_q0: same full graph, up to relabelling
        q0 := Min([x[1] : x in d[4] | x[1] ne p]);
        wp := [y[2] : y in d[4] | y[1] eq p][1]; wq := [y[2] : y in d[4] | y[1] eq q0][1];
        d[4] := [x[1] eq p select <p, [wp[wq[e]] : e in [1..d[2]]]> else x : x in d[4]];
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

// the graph part of a data tuple, without the unit orders
function dg_graph_only(d)
    return <d[1], d[2], d[3], d[4], d[5]>;
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
    // wrong edge termini: noreverse and wrongreversal disconnect the graph; twistedreversal leaves
    // it connected and passes every check on the graph alone, and only the trace of T_p rejects it
    good := DualGraphData(6, 35, 2 : CacheDir := "none");
    for variant in ["noreverse", "wrongreversal", "twistedreversal"] do
        bad := dg_corrupt(good, 2, variant);
        rejected := false;
        try ValidateDualGraphData(6, 35, 2, dg_graph_only(bad)); catch e rejected := true; end try;
        assert rejected eq (variant ne "twistedreversal");
        rejected := false;
        try ValidateDualGraphData(6, 35, 2, bad); catch e rejected := true; end try;
        assert rejected;
    end for;
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
    path := dir cat "/dualgraph_v2_6_35_2";
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
    // a file without the unit orders is rejected on read, even though its graph is right
    Write(path, Sprint(dg_graph_only(d1), "Magma") : Overwrite := true);
    ClearDualGraphCache();
    rejected := false;
    try _ := DualGraphData(6, 35, 2 : CacheDir := dir); catch e rejected := true; end try;
    assert rejected;
    // a file with another format stamp is never read: a v1 file holding a wrong graph is ignored
    System(Sprintf("rm -f '%o'", path));
    Write(dir cat "/dualgraph_v0_6_35_2", "this is not Magma" : Overwrite := true);
    Write(dir cat "/dualgraph_v1_6_35_2", Sprint(dg_graph_only(dg_corrupt(d1, 2, "noreverse")), "Magma") : Overwrite := true);
    ClearDualGraphCache();
    d3 := DualGraphData(6, 35, 2 : CacheDir := dir);
    assert OpenTest(path, "r");
    assert dg_fixcounts(d3[4], d3[2]) eq dg_fixcounts(d1[4], d1[2]);
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
dg_genus_rows := eval Read("tests/_dualgraph_genus_data.txt");

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
    // each deliberately wrong graph (dg_corrupt) breaks some identities.  431 and 397 are regression
    // values measured by this repo at 9a48b57, not independent: all come from W containing a multiple
    // of p, since the genus is blind to the termini otherwise (test_DualGraphBrandtCharpoly checks them)
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
        case variant:
            when "noreverse": assert fails eq 431;
            when "wrongreversal": assert fails eq 397;
            else assert fails gt 0;
        end case;
    end for;
    printf " Done!\n";
end procedure;
test_DualGraphQuotient(dg_genus_rows);

procedure test_DualGraphHyperelliptic()
    printf "Testing canonical model and involution search...";
    // synthetic graphs: rose and banana are hyperelliptic; K4, K33, Petersen, cube are not
    assert IsHyperellipticGraph(1, [<1,1>, <1,1>, <1,1>], 3);
    assert IsHyperellipticGraph(2, [<1,2> : i in [1..4]], 3);
    K4 := [<1,2>, <1,3>, <1,4>, <2,3>, <2,4>, <3,4>];
    assert not IsHyperellipticGraph(4, K4, 3);
    assert not IsHyperellipticGraph(6, [<i,j> : i in [1..3], j in [4..6]], 4);
    petersen := [<i, i mod 5 + 1> : i in [1..5]] cat [<5+i, 5 + ((i+1) mod 5) + 1> : i in [1..5]] cat [<i, i+5> : i in [1..5]];
    petersen := [<Min(e[1], e[2]), Max(e[1], e[2])> : e in petersen];
    t0 := Cputime();
    assert not IsHyperellipticGraph(10, petersen, 6);
    cube := [<a+1, b+1> : a, b in [0..7] | a lt b and #[i : i in [0..2] | ((a div 2^i) mod 2) ne ((b div 2^i) mod 2)] eq 1];
    assert not IsHyperellipticGraph(8, cube, 5);
    assert Cputime(t0) lt 5;
    // canonical model: a theta graph with a pendant path and a subdivided edge collapses to a banana
    nv, edges := CanonicalGraphModel(6, [<1,2>, <1,3>, <3,2>, <1,2>, <2,4>, <4,5>, <5,6>]);
    assert nv eq 2 and #edges eq 3;
    // X_0^6(55)/<w_3, w_110> at p = 3: K4, not hyperelliptic (pilot README)
    hyp, nv, edges := DualGraphTest(6, 55, {1, 3, 110, 330}, 3, 3 : CacheDir := "none");
    assert not hyp;
    assert IsIsomorphic(dg_subdivided(nv, edges), dg_subdivided(4, K4));
    // a wrong genus aborts
    raised := false;
    try _ := DualGraphTest(6, 55, {1, 3, 110, 330}, 4, 3 : CacheDir := "none");
    catch e raised := Position(e`Object, "DUALGRAPH_GENUS") gt 0; end try;
    assert raised;
    // ten known-hyperelliptic curves (pilot hyp.m) pass at every p | D
    for c in [<15,2,{1},3>, <51,2,{1,6},4>, <38,3,{1,3},4>, <57,2,{1,114},3>, <119,1,{1,17},5>,
              <6,23,{1,6},3>, <143,1,{1,13},6>, <146,1,{1,73},3>, <77,2,{1,2,7,14},3>, <33,5,{1,3,11,33},3>] do
        assert GenusShimuraCurveQuotient(c[1], c[2], c[3]) eq c[4];
        for p in PrimeDivisors(c[1]) do
            assert DualGraphTest(c[1], c[2], c[3], c[4], p : CacheDir := "none");
        end for;
    end for;
    printf "Done!\n";
end procedure;
test_DualGraphHyperelliptic();

// a curve with its genus from the independent formula
function dg_mk(D, N, W)
    X := CreateShimuraQuot(D, N, W);
    X`g := GenusShimuraCurveQuotient(D, N, W);
    return X;
end function;

procedure test_FilterByDualGraph()
    printf "Testing FilterByDualGraph...";
    X1 := dg_mk(6, 55, {1, 3, 110, 330});         // ruled at p = 3 only (the graph at p = 2 is hyperelliptic)
    X2 := dg_mk(158, 1, {1, 2});                  // ruled at p = 2 only
    X3 := dg_mk(6, 55, {1, 3, 110, 330}); X3`IsSubhyp := true; X3`IsHyp := true; X3`TestInWhichProved := "hand";
    X4 := CreateShimuraQuot(1, 97, {1}); X4`g := 7;            // D = 1: not applicable
    X5 := CreateShimuraQuot(6, 25, {1}); X5`g := 5;            // N not squarefree: not applicable
    X6 := dg_mk(15, 2, {1});                      // hyperelliptic: no conclusion
    X7 := dg_mk(6, 35, {1, 3, 14, 42});           // ruled at p = 2 and p = 3 (test_DualGraphStankewicz): the tag names 3
    cs := [X1, X2, X3, X4, X5, X6, X7];
    FilterByDualGraph(~cs : CacheDir := "none");
    assert cs[1]`IsSubhyp eq false and cs[1]`IsHyp eq false and cs[1]`TestInWhichProved eq "DualGraph at p = 3";
    assert cs[2]`IsSubhyp eq false and cs[2]`TestInWhichProved eq "DualGraph at p = 2";
    assert cs[3]`IsSubhyp and cs[3]`TestInWhichProved eq "hand";
    assert &and[not assigned cs[i]`IsSubhyp : i in [4, 5, 6]];
    assert cs[7]`IsSubhyp eq false and cs[7]`TestInWhichProved eq "DualGraph at p = 3";
    // a genus inconsistency is not swallowed by the filter
    Y := dg_mk(6, 55, {1, 3, 110, 330}); Y`g := 4;
    ys := [Y];
    raised := false;
    try FilterByDualGraph(~ys : CacheDir := "none");
    catch e raised := Position(e`Object, "DUALGRAPH_GENUS") gt 0; end try;
    assert raised;
    printf "Done!\n";
end procedure;
test_FilterByDualGraph();

procedure test_DualGraphNotSFICertificate()
    printf "Testing that DualGraph verdicts are not SFI certificates...";
    X := CreateShimuraQuot(158, 1, {1, 2}); X`g := GenusShimuraCurveQuotient(158, 1, {1, 2}); X`CurveID := 1;
    cs := [X];
    FilterByDualGraph(~cs : CacheDir := "none");
    assert cs[1]`TestInWhichProved eq "DualGraph at p = 2";
    // ruled by the graph at p = 2 | D; that these primes happen not to certify it is incidental (true
    // today, and a future char-p test could legitimately change it): the textual check below is the guard
    for l in [3, 5, 7] do
        assert not NonHyperellipticAtPrimeCertificate(cs[1], l);
    end for;
    // and the certificate code never consults the dual-graph test
    src := Read("ShimuraQuotients.m");
    a := Position(src, "intrinsic NonHyperellipticAtPrimeCertificate(");
    b := a + Position(src[a..#src], "end intrinsic;");
    assert a gt 0 and Position(src[a..b], "DualGraph") eq 0;
    printf "Done!\n";
end procedure;
test_DualGraphNotSFICertificate();

dg_fiber_raw := eval Read("tests/_dualgraph_fiber_data.txt");

procedure test_DualGraphStankewicz(raw)
    printf "Testing agreement with Stankewicz's fiber...";
    for r in raw do
        D, N, p, ei, al := Explode(r);
        h := Max([x[1] : x in ei]); hh := #ei;
        org := [x[1] : x in ei];
        wp := [x[2] : x in al | x[1] eq p][1];
        assert [x[2] : x in ei] eq [org[wp[e]] : e in [1..hh]];     // terminus = origin of w_p(e)
        alV := [];
        for x in al do
            if x[1] eq p then continue; end if;
            img := [0 : v in [1..h]];
            for e in [1..hh] do img[org[e]] := org[x[2][e]]; end for;
            Append(~alV, <x[1], img>);
        end for;
        theirs := <h, hh, org, al, alV>;
        ValidateDualGraphData(D, N, p, theirs);
        ours := DualGraphData(D, N, p : CacheDir := "none");
        assert dg_same_graphs(ours, theirs, D, N, p);
        // trivial w_p on edges: the comparison rejects it; validation does too (disconnected graph)
        // unless h = 1, where every edge joins the same two vertices and only the quotients differ
        bad := dg_corrupt(theirs, p, "noreverse");
        rejected := false;
        try ValidateDualGraphData(D, N, p, bad); catch e rejected := true; end try;
        assert rejected eq (h gt 1);
        assert not dg_same_graphs(ours, bad, D, N, p);
        // X_0^6(35)/<w_3, w_14>: from their data alone, K4 (not hyperelliptic) at both p = 2 and p = 3
        if <D, N> eq <6, 35> then
            nv, edges := DualGraphQuotientFromData(theirs, p, {1, 3, 14, 42});
            nv, edges := CanonicalGraphModel(nv, edges);
            assert IsIsomorphic(dg_subdivided(nv, edges), dg_subdivided(4, [<1,2>, <1,3>, <1,4>, <2,3>, <2,4>, <3,4>]));
            assert not IsHyperellipticGraph(nv, edges, 3);
        end if;
    end for;
    printf "Done!\n";
end procedure;
test_DualGraphStankewicz(dg_fiber_raw);

// AL permutations from Magma's Brandt module of level M in the algebra of discriminant Dp
function dg_brandt_list(Dp, M)
    B := BrandtModule(Dp, M);
    n := Dimension(B);
    return [<q, [[j : j in [1..n] | A[i][j] ne 0][1] : i in [1..n]]> where A := Matrix(AtkinLehnerOperator(B, q))
            : q in PrimeDivisors(Dp*M)], n;
end function;

procedure test_DualGraphBrandt()
    printf "Testing against Magma's Brandt modules...";
    for c in [<6,35,2>, <10,33,5>, <14,33,7>] do
        D, N, p := Explode(c);
        data := DualGraphData(D, N, p : CacheDir := "none");
        bV, nV := dg_brandt_list(D div p, N);
        bE, nE := dg_brandt_list(D div p, N*p);
        assert nV eq data[1] and nE eq data[2];
        assert dg_fixcounts(bV, nV) eq dg_fixcounts(data[5], nV);
        assert dg_fixcounts(bE, nE) eq dg_fixcounts(data[4], nE);
    end for;
    printf "Done!\n";
end procedure;
test_DualGraphBrandt();

// Independent check of the graph's shape, not only its genus: the unit-weighted +/- adjacency is the
// Brandt matrix T_p of level N in the algebra of discriminant D/p, compared by characteristic
// polynomial (labels differ) with Magma's BrandtModule.  The 13 pilot levels.
procedure test_DualGraphBrandtCharpoly()
    printf "Testing the weighted adjacency against Magma's Brandt matrix T_p...";
    levels := [<6,35,2>, <6,35,3>, <10,33,5>, <14,5,7>, <6,55,3>, <38,3,2>, <38,3,19>, <146,1,2>,
               <119,1,7>, <143,1,11>, <10,21,2>, <22,15,11>, <6,77,2>];
    caught := AssociativeArray();
    for v in ["noreverse", "wrongreversal"] do caught[v] := 0; end for;
    for c in levels do
        D, N, p := Explode(c);
        data := DualGraphData(D, N, p : CacheDir := "none");
        cp := CharacteristicPolynomial(HeckeOperator(BrandtModule(D div p, N), p));
        assert CharacteristicPolynomial(DualGraphBrandtMatrix(data, p)) eq cp;
        for v in ["noreverse", "wrongreversal"] do
            if CharacteristicPolynomial(DualGraphBrandtMatrix(dg_corrupt(data, p, v), p)) ne cp then
                caught[v] +:= 1;
            end if;
        end for;
    end for;
    // wrong termini change T_p at every level with h > 1; at h = 1 it is the 1x1 matrix (p + 1)
    nh := #[c : c in levels | DualGraphData(c[1], c[2], c[3] : CacheDir := "none")[1] gt 1];
    assert nh eq 10 and caught["noreverse"] eq nh and caught["wrongreversal"] eq nh;
    printf "Done!\n";
end procedure;
test_DualGraphBrandtCharpoly();
