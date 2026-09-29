# Dual-Graph Filter Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a pipeline filter that proves X_0^D(N)/W is not hyperelliptic when the dual graph of its special fibre at some p | D has no involution with tree quotient.
**Architecture:** A new package file `DualGraph.m` builds, per (D, N, p), the dual graph of X_0^D(N) at p from left ideal classes of Eichler orders in the definite algebra of discriminant D/p. It caches the graph in memory and on disk. It then quotients the graph by W, reduces it to its canonical model, and searches it for a hyperelliptic involution (unweighted). Two new pipeline stages call it: `FilterByDualGraphStar` in Phase A and `FilterByDualGraph` in Phase B. Both are dealt to workers by level, and Phase B is followed by the closure `UpdateCurvesAfterDualGraph`.
**Tech Stack:** Magma packages via AttachSpec; tests via run_tests.m; bash pipeline scripts.
**Spec:** docs/superpowers/specs/2026-09-29-dual-graph-filter-design.md

## Global Constraints

- Invoke Magma as `magma -b ... < /dev/null > LOG 2>&1`. Always redirect stdin. Never pipe Magma output to `head` or `tail`; read the log file after the run.
- Magma buffers stdout written to a file. An unchanged log does not mean the run has stalled.
- Kill Magma only by PID (`ps -eo pid,etime,command | grep magma.exe`, then `kill <pid>`). Never use `pkill -f magma.exe`.
- Never edit anything under `tests/` while a test suite is running from the same working tree.
- The full suite does not complete on this Mac (it dies at `X0_206_1`). Run only `filename:=` or `target:=` here, and run the whole suite on lava.
- A filter never overwrites a verdict: `if assigned X`IsSubhyp then continue; end if;`.
- A quotient graph whose b_1 differs from `X`g` raises `error "DUALGRAPH_GENUS: ..."`, which stops the stage. It is never caught or skipped.
- `NonHyperellipticAtPrimeCertificate` must never admit a DualGraph verdict. The graph says nothing about reduction at a good prime.
- The test is unweighted only: no unit groups, no stabilizer orders, no edge lengths.
- The filter applies only when D > 1, N is squarefree and gcd(D, N) = 1. Other curves are left untouched and counted as "not applicable".
- `data/dualgraph/` is gitignored and regenerated on lava (user decision, 2026-09-29).
- Implementers commit directly on the branch `dual-graph-filter` (user ruling, 2026-09-29, superseding the earlier "propose commit" wording; treat every "propose commit" step below as "commit"). Never push. Never merge a PR (Eran merges).
- `LeftIdealClasses` labels are not deterministic between calls. Compare two builds only up to isomorphism, never with `eq`, unless both builds come from the same cache file.
- In an eval'd test file, `eval` inside a function or procedure crashes Magma 2.29-4 (segfault, reproduced). Read data files with `eval Read(...)` at the top level of the test file only.

## Review Focus

These five failure modes are implied by the spec but are not exercised by the obvious tests. Each is pinned by a test in the task that owns the code:

1. **Stale disk cache.** A cache file written by older code can be internally consistent and still wrong, and validation cannot detect that. The format stamp `v1` in the file name guards against it. Pinned in Task 2: a file with another stamp (`dualgraph_v0_...`) is never read.
2. **Class-key fallback path.** The full scan in `dg_find` runs only when a key bucket misses, which never happened in any run. So it is dead code in the ordinary tests. Pinned in Task 1, with the full-scan build (`ForceScan := true`) compared to the keyed build by AL fixed-point counts, and in Task 3 up to isomorphism for every W.
3. **Exponential involution search.** Highly symmetric graphs have many involutions and none with tree quotient. Pinned in Task 4: the Petersen graph (g = 6) and the cube graph (g = 5) are both decided in under 5 s in total.
4. **A genus inconsistency swallowed.** A `try` in the filter or the worker could hide `DUALGRAPH_GENUS`. Pinned in Task 5: `FilterByDualGraph` on a curve with a wrong `g` raises `DUALGRAPH_GENUS`.
5. **Corruption that the genus check misses.** Exchanging two origin labels broke only 4 to 15 of 733 genus identities, depending on the random labelling. So validation must catch origin-map errors itself. Pinned in Task 1 and Task 2: a vertex action that the origin map does not intertwine is rejected by `ValidateDualGraphData`, both on build and when read from disk.

## Task 0: Worktree and branch

**Files:** none (git only).
**Interfaces:** none.

The pipeline files this plan modifies (`run_pipeline.sh`, `run_parallel_filter.sh`, `parallel_filter_worker.m` with its level grouping, `run_sequential_stage.m`, `docs/RUNNING_PIPELINE.md`, and `NonHyperellipticAtPrimeCertificate`) exist only on `integration`. Checked on 2026-09-29:

* `git log --oneline -3 integration` gives `00b1f0b Merge fix/sfi-char-p into integration (#51 review follow-ups)`, `ba9a09f`, `dc353d6`.
* `main` (`989a130`) is an ancestor of `integration`, and `integration` is 481 commits ahead.
* `main` has no `Twisted` stage in `run_pipeline.sh`, and no `NonHyperellipticAtPrimeCertificate`.

So branch from `integration`. CLAUDE.md says "two branches, that is all". Ask the user before creating the branch.

- [ ] **Step 1: Confirm the base.** Run `git -C /Users/sachihashimoto/Repos/ShimuraCurveALQuotients log --oneline -3 integration`. Expected: first line `00b1f0b Merge fix/sfi-char-p into integration ...`. If it differs, report the new head before continuing.
- [ ] **Step 2: Ask the user** to approve `git worktree add worktrees/dual-graph -b dual-graph-filter integration`, run from the repo root. Once approved, run it. Every later path is relative to `worktrees/dual-graph/`, and every command runs from there.
- [ ] **Step 3: Check the environment.** Run `magma -b filename:=tests/Kappa0.m run_tests.m < /dev/null > /tmp/dg_task0.log 2>&1`. Expected in the log: `Kappa0.m: Testing Kappa0...Done!` and `Success!`.

## Task 1: Data layer (`DualGraph.m`) with in-memory store and validation

**Files:**
- Create: `DualGraph.m`
- Modify: `Caching.m:8`. After `sfi_certificates := NewStore();`, add the store.
- Modify: `ShimuraQuotients.spec:31`. After `special_fiber_cm.m`, add `DualGraph.m`.
- Create: `tests/DualGraph.m`

**Interfaces:**
- Consumes: `GenusShimuraCurveQuotient(D::RngIntElt, N::RngIntElt, als::SetEnum) -> RngIntElt` (`ShimuraQuotients.m:453`), and `GetCache`, `SetCache`, `NewStore` (`Caching.m`).
- Produces:
  - `DualGraphApplicable(D::RngIntElt, N::RngIntElt) -> BoolElt`
  - `DualGraphData(D::RngIntElt, N::RngIntElt, p::RngIntElt : CacheDir := "", ForceScan := false) -> Tup`, returning `<h, h', org, alE, alV>`. Here alE and alV are sequences of `<q, images>`.
  - `ValidateDualGraphData(D::RngIntElt, N::RngIntElt, p::RngIntElt, data::Tup)`
  - `ClearDualGraphCache()`

- [ ] **Step 1: Write the failing test.** Create `tests/DualGraph.m` with this content:

<!-- file: tests/DualGraph.m create -->
```magma
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
```

- [ ] **Step 2: Run it and watch it fail.** Run `magma -b filename:=tests/DualGraph.m run_tests.m < /dev/null > /tmp/dg_task1.log 2>&1`. Expected in the log: `DualGraph.m: Fail!`, with `Identifier 'DualGraphData' has not been declared or assigned`.

- [ ] **Step 3: Add the store to `Caching.m`.** After line 8 (`sfi_certificates := NewStore();`), insert:

```magma
// DualGraphData fibers of X_0^D(N) at p | D, keyed by <D, N, p>
dual_graphs := NewStore();
```

- [ ] **Step 4: Register the file.** In `ShimuraQuotients.spec`, after the line `special_fiber_cm.m` (line 31), add the line `DualGraph.m`.

- [ ] **Step 5: Create `DualGraph.m`.**

<!-- file: DualGraph.m create -->
```magma
// Dual graph of the special fibre at p | D of X_0^D(N)/W, and the graph-hyperellipticity test.
//
// Y = X_0^D(N)/W, D > 1, N squarefree, gcd(D, N) = 1, p | D.  Y over Z_p is a Mumford curve whose
// stable special fibre has dual graph the minimised quotient graph G_W (Padurariu--Saia,
// arXiv:2509.25368, Thm 2.6, after Ogg 1985 and Bertolini--Darmon 1996).  gon(graph) <= gon(Y)
// (Baker 2008).  A 2-edge-connected graph of genus >= 2 has gonality 2 iff it has an involution
// whose quotient is a tree (Baker--Norine 2009, Chan 2013), i.e. with chi(Fix) = g + 1.  Lengths
// are ignored (unweighted): an isometric involution is in particular a combinatorial one, so "no
// combinatorial involution" is sound for every choice of lengths.
//
// Graph G (Kurihara; Padurariu--Saia 2.3).  B' = definite algebra of discriminant D/p, O an Eichler
// order of level N, OO the level-Np suborder.  Vertices (v, +) and (v, -) for the h left ideal
// classes v of O; one edge e for each of the h' left ideal classes of OO, joining (org[e], +) to
// (org[w_p e], -), where org[e] is the class of O*I_e.  w_q (q | D'N) acts on classes of both
// levels by I -> P_q I (P_q the two-sided prime of norm q) and preserves the copies; w_p acts on
// edges by P_p and swaps the copies.  Specification: docs/superpowers/specs/2026-09-29-dual-graph-filter-design.md.

import "Caching.m" : dual_graphs, SetCache, GetCache;

// ---------- ideal-class lookup ----------

// Class invariant of a left ideal: canonical Gram matrix of its norm form scaled by 1/Norm(I),
// unchanged by I -> I*alpha.  Only a bucket key: every lookup is confirmed by IsIsomorphic.
function dg_key(I)
    return MinkowskiGramReduction(GramMatrix(I)/Norm(I) : Canonical := true);
end function;

function dg_buckets(reps)
    b := AssociativeArray();
    for i->I in reps do
        k := dg_key(I);
        if IsDefined(b, k) then Append(~b[k], i); else b[k] := [i]; end if;
    end for;
    return b;
end function;

// index of the class of J among reps; falls back to a full scan if the key bucket has no match
function dg_find(J, reps, buckets)
    ok, cand := IsDefined(buckets, dg_key(J));
    if ok then
        hits := [i : i in cand | IsIsomorphic(J, reps[i])];
        if #hits eq 1 then return hits[1]; end if;
    end if;
    vprintf ShimuraQuotients, 1 : "DualGraph: class key missed, scanning all %o classes\n", #reps;
    hits := [i : i in [1..#reps] | IsIsomorphic(J, reps[i])];
    assert #hits eq 1;
    return hits[1];
end function;

function dg_prime(O, q)
    P := PrimeIdeal(O, q);
    assert Norm(P) eq q and IsTwoSidedIdeal(P);
    return P;
end function;

// image sequence of w_q in a list of pairs <q, images>
function dg_img(list, q)
    i := [j : j in [1..#list] | list[j][1] eq q];
    assert #i eq 1;
    return list[i[1]][2];
end function;

// ---------- the data layer ----------

function dg_build(D, N, p, scan)
    Dp := D div p;
    Q := QuaternionAlgebra(Dp);
    O := QuaternionOrder(Q, N);
    OO := Order(O, p);
    L := LeftIdealClasses(O);
    LL := LeftIdealClasses(OO);
    // scan: empty buckets, so every lookup takes the full-scan path (used by tests of that path)
    bL := scan select AssociativeArray() else dg_buckets(L);
    bLL := scan select AssociativeArray() else dg_buckets(LL);
    org := [dg_find(lideal< O | Basis(II) >, L, bL) : II in LL];
    alE := [];
    alV := [];
    for q in PrimeDivisors(Dp*N*p) do
        P := dg_prime(OO, q);
        Append(~alE, <q, [dg_find(P*II, LL, bLL) : II in LL]>);
        if q ne p then
            P0 := dg_prime(O, q);
            Append(~alV, <q, [dg_find(P0*I, L, bL) : I in L]>);
        end if;
    end for;
    return <#L, #LL, org, alE, alV>;
end function;

intrinsic ValidateDualGraphData(D::RngIntElt, N::RngIntElt, p::RngIntElt, data::Tup)
{Raises an error unless data is a consistent dual graph for X_0^D(N) at p: the Atkin-Lehner maps
are involutions of the right sets, the origin map is equivariant for every w_q with q | DN/p, and
h' - 2h + 1 is the genus of X_0^D(N).  Run on every build and on every disk-cache read.}
    require #data eq 5 : "data must be <h, h', org, alE, alV>";
    h := data[1]; hh := data[2]; org := data[3];
    assert #org eq hh and &and[v in [1..h] : v in org];
    assert {x[1] : x in data[4]} eq Set(PrimeDivisors(D*N));
    assert {x[1] : x in data[5]} eq Set(PrimeDivisors((D div p)*N));
    for x in data[4] do
        img := x[2];
        assert Sort(img) eq [1..hh] and &and[img[img[e]] eq e : e in [1..hh]];
    end for;
    for x in data[5] do
        q := x[1]; img := x[2]; eimg := dg_img(data[4], q);
        assert Sort(img) eq [1..h] and &and[img[img[v]] eq v : v in [1..h]];
        assert &and[org[eimg[e]] eq img[org[e]] : e in [1..hh]];       // origin map is equivariant
    end for;
    // b_1 of the full graph (2h vertices, h' edges, connected) is the genus of X_0^D(N)
    assert hh - 2*h + 1 eq GenusShimuraCurveQuotient(D, N, {Integers() | 1});
end intrinsic;

intrinsic ClearDualGraphCache()
{Empties the in-memory DualGraphData cache (the disk cache is untouched).}
    StoreClear(dual_graphs);
end intrinsic;

intrinsic DualGraphApplicable(D::RngIntElt, N::RngIntElt) -> BoolElt
{True iff the dual-graph test applies to X_0^D(N)/W: D > 1, N squarefree, gcd(D, N) = 1.}
    return D gt 1 and IsSquarefree(N) and GCD(D, N) eq 1;
end intrinsic;

intrinsic DualGraphData(D::RngIntElt, N::RngIntElt, p::RngIntElt : CacheDir := "", ForceScan := false) -> Tup
{The dual graph of the special fibre at p of X_0^D(N), as <h, h', org, alE, alV>: h left ideal
classes of an Eichler order of level N in the definite algebra of discriminant D/p (one vertex in
each of two copies), h' classes of level Np (the edges; edge e joins (org[e], +) to
(org[alE_p[e]], -)), and the Atkin-Lehner involutions as pairs <q, images>, on edges for q | DN and
on vertices for q | DN/p.  Cached per <D, N, p> for the session.  CacheDir names the disk cache,
which this version does not read or write yet.  ForceScan builds afresh (no caches) identifying
every class by a full scan instead of the key buckets.}
    require DualGraphApplicable(D, N) : "needs D > 1, N squarefree, gcd(D, N) = 1";
    require IsPrime(p) and D mod p eq 0 : "p must be a prime dividing D";
    if ForceScan then
        data := dg_build(D, N, p, true);
        ValidateDualGraphData(D, N, p, data);
        return data;
    end if;
    key := <D, N, p>;
    b, data := GetCache(key, dual_graphs);
    if b then return data; end if;
    data := dg_build(D, N, p, false);
    ValidateDualGraphData(D, N, p, data);
    SetCache(key, data, dual_graphs);
    return data;
end intrinsic;
```

- [ ] **Step 6: Run it again.** Run the same command as in Step 2. Expected: `Testing DualGraphData...Done!` and `DualGraph.m: Success!`. This takes about 5 s.
- [ ] **Step 7: Propose a commit** of `DualGraph.m`, `Caching.m`, `ShimuraQuotients.spec` and `tests/DualGraph.m`, with the message:
  `DualGraph.m: dual graph of X_0^D(N) at p | D from ideal classes, validated on build`

## Task 2: Disk cache (format-stamped, tmp + mv, re-validated on read)

**Files:**
- Modify: `DualGraph.m`. Insert the disk-cache functions before `intrinsic ClearDualGraphCache`, and replace the whole `intrinsic DualGraphData` from Task 1.
- Modify: `.gitignore`. Append at the end (line 28).
- Test: `tests/DualGraph.m` (append).

**Interfaces:**
- Consumes: `ValidateDualGraphData` (Task 1).
- Produces: the same `DualGraphData` signature. The disk cache is now live: `CacheDir` defaults to `$DUALGRAPH_CACHE_DIR`, else `data/dualgraph`, and the value `"none"` disables it. The file is `<dir>/dualgraph_v1_<D>_<N>_<p>`.

- [ ] **Step 1: Write the failing test.** Append to `tests/DualGraph.m`:

<!-- file: tests/DualGraph.m append -->
```magma
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
```

- [ ] **Step 2: Run it and watch it fail.** Run `magma -b filename:=tests/DualGraph.m run_tests.m < /dev/null > /tmp/dg_task2.log 2>&1`. Expected: `Testing DualGraph disk cache...` followed by `DualGraph.m: Fail!` with `Assertion failed`, because no file has been written.

- [ ] **Step 3: Add the disk-cache functions.** In `DualGraph.m`, insert this block before `intrinsic ClearDualGraphCache()`:

<!-- file: DualGraph.m insert-before "intrinsic ClearDualGraphCache()" -->
```magma
// ---------- disk cache: <dir>/dualgraph_<format>_<D>_<N>_<p> ----------

function dg_cache_dir(CacheDir)
    if CacheDir ne "" then return CacheDir; end if;
    d := GetEnv("DUALGRAPH_CACHE_DIR");
    return d eq "" select "data/dualgraph" else d;
end function;

// The format stamp is part of the name: bump it whenever dg_build or the tuple layout changes, so a
// cache written by older code is never read (validation cannot catch consistent-but-stale data).
DG_FORMAT := "v1";

function dg_cache_path(dir, D, N, p)
    return Sprintf("%o/dualgraph_%o_%o_%o_%o", dir, DG_FORMAT, D, N, p);
end function;

// Write to a per-process temporary name, then rename: a killed worker never leaves a partial file.
procedure dg_cache_write(dir, D, N, p, data)
    System(Sprintf("mkdir -p '%o'", dir));
    path := dg_cache_path(dir, D, N, p);
    tmp := Sprintf("%o.tmp.%o", path, Getpid());
    Write(tmp, Sprint(data, "Magma") : Overwrite := true);
    System(Sprintf("mv '%o' '%o'", tmp, path));
end procedure;

```

- [ ] **Step 4: Replace `DualGraphData`.** Replace the whole `intrinsic DualGraphData ... end intrinsic;` with:

<!-- file: DualGraph.m replace "intrinsic DualGraphData" -->
```magma
intrinsic DualGraphData(D::RngIntElt, N::RngIntElt, p::RngIntElt : CacheDir := "", ForceScan := false) -> Tup
{The dual graph of the special fibre at p of X_0^D(N), as <h, h', org, alE, alV>: h left ideal
classes of an Eichler order of level N in the definite algebra of discriminant D/p (one vertex in
each of two copies), h' classes of level Np (the edges; edge e joins (org[e], +) to
(org[alE_p[e]], -)), and the Atkin-Lehner involutions as pairs <q, images>, on edges for q | DN and
on vertices for q | DN/p.  Cached per <D, N, p> for the session and on disk in CacheDir (default
$DUALGRAPH_CACHE_DIR, else data/dualgraph; "none" disables the disk cache).  A disk entry is
re-validated when read, so a corrupt file raises an error instead of giving a verdict.  ForceScan
builds afresh (no caches) identifying every class by a full scan instead of the key buckets.}
    require DualGraphApplicable(D, N) : "needs D > 1, N squarefree, gcd(D, N) = 1";
    require IsPrime(p) and D mod p eq 0 : "p must be a prime dividing D";
    if ForceScan then
        data := dg_build(D, N, p, true);
        ValidateDualGraphData(D, N, p, data);
        return data;
    end if;
    key := <D, N, p>;
    b, data := GetCache(key, dual_graphs);
    if b then return data; end if;
    dir := dg_cache_dir(CacheDir);
    path := dg_cache_path(dir, D, N, p);
    if dir ne "none" and OpenTest(path, "r") then
        data := eval Read(path);
        ValidateDualGraphData(D, N, p, data);
    else
        data := dg_build(D, N, p, false);
        ValidateDualGraphData(D, N, p, data);
        if dir ne "none" then dg_cache_write(dir, D, N, p, data); end if;
    end if;
    SetCache(key, data, dual_graphs);
    return data;
end intrinsic;
```

- [ ] **Step 5: Gitignore the cache.** Append to `.gitignore`:

```
# DualGraphData disk cache (DualGraph.m); regenerated on lava, never committed
data/dualgraph/
```

- [ ] **Step 6: Run it again.** Run the command from Step 2. Expected: both `Done!` lines, and `DualGraph.m: Success!`.
- [ ] **Step 7: Propose a commit** with the message:
  `DualGraph.m: format-stamped disk cache, written by tmp+mv, re-validated on read`

## Task 3: Quotient graph for W and the genus identities

**Files:**
- Modify: `DualGraph.m`. Append the quotient intrinsic at the end.
- Create: `tests/dualgraph_genus_data.txt`
- Test: `tests/DualGraph.m` (append).

**Interfaces:**
- Consumes: `DualGraphData` (Task 2).
- Produces: `DualGraphQuotientFromData(data::Tup, p::RngIntElt, W::SetEnum) -> RngIntElt, SeqEnum`. The results are the number of vertex orbits that carry an edge, and the edges `<u, v>` with u ≤ v. Half-edges are removed and leaves are kept.

- [ ] **Step 1: Generate the data file.** Run from the worktree root:

```bash
python3 - <<'EOF'
import re
from math import gcd
src='/Users/sachihashimoto/Repos/ShimuraCurveALQuotients/handoff_2026-09-29/graph-pilot/genustest.m'
recs=re.findall(r'\[\* (\d+), (\d+), (\d+), \[([^\]]*)\], (\d+) \*\]', open(src).read())
def sqf(n): return all(n % (q*q) for q in range(2, int(n**0.5)+1))
rows=[f"<{D}, {N}, {{ {W} }}, {g}>" for _,D,N,W,g in recs if sqf(int(N)) and gcd(int(D),int(N))==1]
open('tests/dualgraph_genus_data.txt','w').write("// (D, N, W, g) from handoff_2026-09-29/graph-pilot/genustest.m, N squarefree; the genus g is\n// GenusShimuraCurveQuotient(D, N, W).  Every (D, N, W, p), p | D, is a genus identity b_1(G_W) = g.\n[\n"+",\n".join(rows)+"\n]\n")
print(len(rows))
EOF
```

  Expected output: `412`.

- [ ] **Step 2: Write the failing test.** Append to `tests/DualGraph.m`:

<!-- file: tests/DualGraph.m append -->
```magma
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
```

- [ ] **Step 3: Run it and watch it fail.** Run `magma -b filename:=tests/DualGraph.m run_tests.m < /dev/null > /tmp/dg_task3.log 2>&1`. Expected: `DualGraph.m: Fail!`, with `Identifier 'DualGraphQuotientFromData' has not been declared or assigned`.

- [ ] **Step 4: Append to `DualGraph.m`.**

<!-- file: DualGraph.m append -->
```magma

// ---------- quotient graph for W ----------

intrinsic DualGraphQuotientFromData(data::Tup, p::RngIntElt, W::SetEnum) -> RngIntElt, SeqEnum
{The quotient of the dual graph `data` (from DualGraphData at p) by the Atkin-Lehner group W (a set
of Hall divisors of DN), with half-edges removed and leaves kept.  Returns the number of vertex
orbits that carry an edge and the edges as pairs <u, v>, u <= v.}
    h := data[1]; hh := data[2]; org := data[3];
    wp := dg_img(data[4], p);
    // vertex x in [1..2h]: label ((x-1) mod h) + 1, copy + if x <= h, - otherwise
    function vact(m, x)
        v := ((x-1) mod h) + 1; s := (x-1) div h;
        for q in PrimeDivisors(m) do
            if q eq p then s := 1 - s; else v := dg_img(data[5], q)[v]; end if;
        end for;
        return v + s*h;
    end function;
    function eact(m, e)
        for q in PrimeDivisors(m) do e := dg_img(data[4], q)[e]; end for;
        return e;
    end function;
    vorb := [0 : x in [1..2*h]]; nv := 0;
    for x in [1..2*h] do
        if vorb[x] ne 0 then continue; end if;
        nv +:= 1;
        for m in W do vorb[vact(m, x)] := nv; end for;
    end for;
    edges := [];
    seen := {};
    for e in [1..hh] do
        if e in seen then continue; end if;
        seen join:= {eact(m, e) : m in W};
        // an m with p | m that fixes e swaps its two ends: the orbit is a half-edge, dropped
        if exists{m : m in W | m mod p eq 0 and eact(m, e) eq e} then continue; end if;
        u := vorb[org[e]]; v := vorb[org[wp[e]] + h];
        Append(~edges, <Min(u, v), Max(u, v)>);
    end for;
    used := #edges eq 0 select 0 else #({x[1] : x in edges} join {x[2] : x in edges});
    return used, edges;
end intrinsic;
```

- [ ] **Step 5: Run it again.** Run the command from Step 3. Expected: `Testing DualGraphQuotientFromData... noreverse:431 wrongreversal:397 trivialvertexAL:N Done!` and `DualGraph.m: Success!`, in about 15 s. The counts 431 and 397 were measured on 2026-09-29 and are the same in every session. N varies between runs (279–305 observed), because that corruption breaks equivariance and so the choice of orbit representative matters. Only `fails gt 0` is asserted.
- [ ] **Step 6: Propose a commit** with the message:
  `DualGraph.m: quotient graph by W (half-edges dropped); 822 genus identities and three negative controls`

## Task 4: Canonical model, involution search and `DualGraphTest`

**Files:**
- Modify: `DualGraph.m` (append).
- Test: `tests/DualGraph.m` (append).

**Interfaces:**
- Consumes: `DualGraphData` and `DualGraphQuotientFromData`.
- Produces:
  - `CanonicalGraphModel(nv::RngIntElt, edges::SeqEnum) -> RngIntElt, SeqEnum`
  - `IsHyperellipticGraph(nv::RngIntElt, edges::SeqEnum, g::RngIntElt) -> BoolElt, SeqEnum`
  - `DualGraphTest(D::RngIntElt, N::RngIntElt, W::SetEnum, g::RngIntElt, p::RngIntElt : CacheDir := "") -> BoolElt, RngIntElt, SeqEnum`

- [ ] **Step 1: Write the failing test.** Append to `tests/DualGraph.m`:

<!-- file: tests/DualGraph.m append -->
```magma
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
```

- [ ] **Step 2: Run it and watch it fail.** Run `magma -b filename:=tests/DualGraph.m run_tests.m < /dev/null > /tmp/dg_task4.log 2>&1`. Expected: `Fail!`, with `Identifier 'IsHyperellipticGraph' has not been declared or assigned`.

- [ ] **Step 3: Append to `DualGraph.m`.** This is a port of `handoff_2026-09-29/graph-pilot/graphtest.m:17–140`, with lengths dropped (unweighted), so multiplicities are integers:

<!-- file: DualGraph.m append -->
```magma

// ---------- canonical model and involution search ----------

function dg_relabel(edges)
    usedv := Sort(Setseq({e[1] : e in edges} join {e[2] : e in edges}));
    idx := AssociativeArray();
    for i->v in usedv do idx[v] := i; end for;
    return #usedv, [<idx[e[1]], idx[e[2]]> : e in edges];
end function;

function dg_connected_without(nv, edges, skip)
    seen := {1}; frontier := [1];
    while #frontier gt 0 do
        x := frontier[#frontier]; Prune(~frontier);
        for i->e in edges do
            if i eq skip then continue; end if;
            if e[1] eq x and not (e[2] in seen) then Include(~seen, e[2]); Append(~frontier, e[2]); end if;
            if e[2] eq x and not (e[1] in seen) then Include(~seen, e[1]); Append(~frontier, e[1]); end if;
        end for;
    end while;
    return #seen eq nv;
end function;

intrinsic CanonicalGraphModel(nv::RngIntElt, edges::SeqEnum) -> RngIntElt, SeqEnum
{Canonical model of a connected multigraph given by edges <u, v> (loops allowed): repeatedly delete
leaves, contract bridges and merge the two edges at a valence-2 vertex.  The genus is unchanged.}
    if #edges eq 0 then return 0, edges; end if;
    g0 := #edges - nv + 1;
    changed := true;
    while changed do
        changed := false;
        nv, edges := dg_relabel(edges);
        if #edges eq 0 then break; end if;
        for v in [1..nv] do                                   // leaves
            inc := [i : i->e in edges | e[1] eq v or e[2] eq v];
            if #inc eq 1 and edges[inc[1]][1] ne edges[inc[1]][2] then
                Remove(~edges, inc[1]); changed := true; break;
            end if;
        end for;
        if changed then continue; end if;
        for i->e in edges do                                  // bridges: contract
            if e[1] ne e[2] and not dg_connected_without(nv, edges, i) then
                u := e[1]; v := e[2];
                Remove(~edges, i);
                edges := [<f[1] eq v select u else f[1], f[2] eq v select u else f[2]> : f in edges];
                edges := [<Min(f[1], f[2]), Max(f[1], f[2])> : f in edges];
                changed := true; break;
            end if;
        end for;
        if changed then continue; end if;
        for v in [1..nv] do                                   // valence-2 vertices: merge
            inc := [i : i->e in edges | e[1] eq v or e[2] eq v];
            if #inc eq 2 and forall{i : i in inc | edges[i][1] ne edges[i][2]} then
                e1 := edges[inc[1]]; e2 := edges[inc[2]];
                a := e1[1] eq v select e1[2] else e1[1];
                b := e2[1] eq v select e2[2] else e2[1];
                edges := [edges[i] : i in [1..#edges] | not (i in inc)] cat [<Min(a, b), Max(a, b)>];
                changed := true; break;
            end if;
        end for;
    end while;
    if #edges eq 0 then return 0, edges; end if;
    nv, edges := dg_relabel(edges);
    assert #edges - nv + 1 eq g0;
    return nv, edges;
end intrinsic;

// Backtracking over involutions pi of the vertices preserving edge multiplicities.  At a full
// assignment, the largest chi(Fix) over the compatible edge maps is: fixed vertices, + loops at
// fixed vertices (each reversed, midpoint fixed), - (mult mod 2) for each pair of distinct fixed
// vertices (parallel edges swapped in pairs, one left fixed), + mult for each swapped pair a <-> b
// (each such edge reversed).  chi(Fix) <= g + 1 always, with equality iff the quotient is a tree.
function dg_search(n, pi, mult, g)
    u := 0;
    for i in [1..n] do if pi[i] eq 0 then u := i; break; end if; end for;
    if u eq 0 then
        chi := #[i : i in [1..n] | pi[i] eq i];
        for a in [1..n] do
            for b in [a..n] do
                m := mult(a, b);
                if a eq b then
                    if pi[a] eq a then chi +:= m; end if;
                elif pi[a] eq a and pi[b] eq b then
                    chi -:= m mod 2;
                elif pi[a] eq b then
                    chi +:= m;
                end if;
            end for;
        end for;
        assert chi le g + 1;
        return chi eq g + 1, pi;
    end if;
    for v in [u] cat [w : w in [u+1..n] | pi[w] eq 0] do
        pi2 := pi; pi2[u] := v; pi2[v] := u;
        asg := [i : i in [1..n] | pi2[i] ne 0];
        if forall{<a, b> : a, b in asg | b lt a or mult(a, b) eq mult(pi2[a], pi2[b])} then
            found, w := dg_search(n, pi2, mult, g);
            if found then return true, w; end if;
        end if;
    end for;
    return false, pi;
end function;

intrinsic IsHyperellipticGraph(nv::RngIntElt, edges::SeqEnum, g::RngIntElt) -> BoolElt, SeqEnum
{For a 2-edge-connected multigraph of genus g >= 2 (edges <u, v> on [1..nv]), true iff it has an
involution whose quotient is a tree; if so, also returns it as the vertex images.}
    require g ge 2 and #edges - nv + 1 eq g : "needs a graph of genus g >= 2";
    M := AssociativeArray();
    for e in edges do
        k := <Min(e[1], e[2]), Max(e[1], e[2])>;
        M[k] := IsDefined(M, k) select M[k] + 1 else 1;
    end for;
    mult := func< a, b | IsDefined(M, <Min(a, b), Max(a, b)>) select M[<Min(a, b), Max(a, b)>] else 0 >;
    return dg_search(nv, [0 : i in [1..nv]], mult, g);
end intrinsic;

intrinsic DualGraphTest(D::RngIntElt, N::RngIntElt, W::SetEnum, g::RngIntElt, p::RngIntElt : CacheDir := "") -> BoolElt, RngIntElt, SeqEnum
{For Y = X_0^D(N)/W of genus g >= 2 and p | D: false if the dual graph of Y at p has no involution
with tree quotient, which proves Y is not hyperelliptic; true otherwise (no conclusion).  Also
returns the canonical model (vertex count, edges).  Raises DUALGRAPH_GENUS if the quotient graph
does not have genus g.}
    require g ge 2 : "needs g >= 2";
    data := DualGraphData(D, N, p : CacheDir := CacheDir);
    nv, edges := DualGraphQuotientFromData(data, p, W);
    b1 := #edges eq 0 select 0 else #edges - nv + 1;
    if b1 ne g then
        error Sprintf("DUALGRAPH_GENUS: X_0^%o(%o)/%o at p = %o has graph genus %o, curve genus %o", D, N, W, p, b1, g);
    end if;
    nv, edges := CanonicalGraphModel(nv, edges);
    hyp := IsHyperellipticGraph(nv, edges, g);
    return hyp, nv, edges;
end intrinsic;
```

- [ ] **Step 4: Run it again.** Expected: `Testing canonical model and involution search...Done!` and `Success!`.
- [ ] **Step 5: Propose a commit** with the message:
  `DualGraph.m: canonical model and tree-quotient involution search (unweighted); K4 and controls`

## Task 5: `FilterByDualGraph` (verdicts, odd primes first, not-applicable count)

**Files:**
- Modify: `DualGraph.m` (append).
- Test: `tests/DualGraph.m` (append).

**Interfaces:**
- Consumes: `DualGraphTest` and `DualGraphApplicable`.
- Produces:
  - `DualGraphVerdict(X::ShimuraQuot : CacheDir := "") -> BoolElt, RngIntElt`
  - `FilterByDualGraph(~curves::SeqEnum : CacheDir := "")`, which prints one line:
    `DualGraph summary: tested T, ruled R (witness p = 2: R2), not applicable A, already decided C, genus < 3 L`.

- [ ] **Step 1: Write the failing test.** Append to `tests/DualGraph.m`:

<!-- file: tests/DualGraph.m append -->
```magma
// a curve with its genus from the independent formula
function dg_mk(D, N, W)
    X := CreateShimuraQuot(D, N, W);
    X`g := GenusShimuraCurveQuotient(D, N, W);
    return X;
end function;

procedure test_FilterByDualGraph()
    printf "Testing FilterByDualGraph...";
    X1 := dg_mk(6, 55, {1, 3, 110, 330});         // ruled at p = 2 and p = 3: the tag names 3
    X2 := dg_mk(158, 1, {1, 2});                  // ruled at p = 2 only
    X3 := dg_mk(6, 55, {1, 3, 110, 330}); X3`IsSubhyp := true; X3`IsHyp := true; X3`TestInWhichProved := "hand";
    X4 := CreateShimuraQuot(1, 97, {1}); X4`g := 7;            // D = 1: not applicable
    X5 := CreateShimuraQuot(6, 25, {1}); X5`g := 5;            // N not squarefree: not applicable
    X6 := dg_mk(15, 2, {1});                      // hyperelliptic: no conclusion
    cs := [X1, X2, X3, X4, X5, X6];
    FilterByDualGraph(~cs : CacheDir := "none");
    assert cs[1]`IsSubhyp eq false and cs[1]`IsHyp eq false and cs[1]`TestInWhichProved eq "DualGraph at p = 3";
    assert cs[2]`IsSubhyp eq false and cs[2]`TestInWhichProved eq "DualGraph at p = 2";
    assert cs[3]`IsSubhyp and cs[3]`TestInWhichProved eq "hand";
    assert &and[not assigned cs[i]`IsSubhyp : i in [4, 5, 6]];
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
```

- [ ] **Step 2: Run it and watch it fail.** Run `magma -b filename:=tests/DualGraph.m run_tests.m < /dev/null > /tmp/dg_task5.log 2>&1`. Expected: `Fail!`, with `Identifier 'FilterByDualGraph' has not been declared or assigned`.

- [ ] **Step 3: Append to `DualGraph.m`.** It follows the verdict conventions of `FilterByTrace` (`ShimuraQuotients.m:773–787`).

<!-- file: DualGraph.m append -->
```magma

// ---------- the filter ----------

// odd primes first, ascending, then 2: a witness p = 2 then means no odd prime of D sufficed
function dg_primes(D)
    ps := PrimeDivisors(D);
    return [q : q in ps | q ne 2] cat [q : q in ps | q eq 2];
end function;

intrinsic DualGraphVerdict(X::ShimuraQuot : CacheDir := "") -> BoolElt, RngIntElt
{False and a prime p | D if the dual graph of X at p proves X is not hyperelliptic; true if no
prime of D does (no conclusion).  Needs X`g >= 3 and DualGraphApplicable(X`D, X`N).}
    require X`g ge 3 : "needs genus >= 3";
    require DualGraphApplicable(X`D, X`N) : "needs D > 1, N squarefree, gcd(D, N) = 1";
    for p in dg_primes(X`D) do
        hyp := DualGraphTest(X`D, X`N, X`W, X`g, p : CacheDir := CacheDir);
        if not hyp then return false, p; end if;
    end for;
    return true, _;
end intrinsic;

intrinsic FilterByDualGraph(~curves::SeqEnum : CacheDir := "")
{Marks as not hyperelliptic every undecided curve of genus >= 3 whose dual graph at some p | D has
no involution with tree quotient.  Curves with D = 1 or N not squarefree are not applicable and are
left unchanged; decided curves are never overwritten.  Prints one DualGraph summary line.}
    tested := 0; ruled := 0; ruled2 := 0; notapp := 0; decided := 0; lowg := 0;
    for i->X in curves do
        if assigned X`IsSubhyp then decided +:= 1; continue; end if;
        if X`g lt 3 then lowg +:= 1; continue; end if;
        if not DualGraphApplicable(X`D, X`N) then notapp +:= 1; continue; end if;
        tested +:= 1;
        hyp, p := DualGraphVerdict(X : CacheDir := CacheDir);
        if not hyp then
            curves[i]`IsSubhyp := false;
            curves[i]`IsHyp := false;
            curves[i]`TestInWhichProved := Sprintf("DualGraph at p = %o", p);
            ruled +:= 1;
            if p eq 2 then ruled2 +:= 1; end if;
        end if;
    end for;
    printf "DualGraph summary: tested %o, ruled %o (witness p = 2: %o), not applicable %o, already decided %o, genus < 3 %o\n",
        tested, ruled, ruled2, notapp, decided, lowg;
end intrinsic;
```

- [ ] **Step 4: Run it again.** Expected in the log: `Testing FilterByDualGraph...DualGraph summary: tested 3, ruled 2 (witness p = 2: 1), not applicable 2, already decided 1, genus < 3 0`, then `Done!` and `Success!`.
- [ ] **Step 5: Propose a commit** with the message:
  `FilterByDualGraph: verdicts tagged "DualGraph at p = <p>", odd primes first, never overwrites`

## Task 6: Pipeline wiring

**Files:**
- Modify: `workingcode.m:33`. After `<"HHProposition1", CheckHHProposition1>,` insert `<"FilterByDualGraphStar", FilterByDualGraph>,`.
- Modify: `workingcode.m:57`. After `<"UpdateCurves2", UpdateCurves>,` insert two lines, `<"FilterByDualGraph", FilterByDualGraph>,` and then `<"UpdateCurvesAfterDualGraph", UpdateCurves>,`.
- Modify: `run_pipeline.sh`, at lines 19–36 (header), 126 (Phase A) and 161–164 (Phase B).
- Modify: `run_parallel_filter.sh:32–105`.
- Modify: `parallel_filter_worker.m:55–58`, `:80–81` and `:111`.
- Modify: `run_sequential_stage.m:62–66`.
- Modify: `Utils.m:25` (`CurveCostProxy`).
- Modify: `reconstruct_attribution.m:34` and `:48`, `analysis_stages.m:20` and `:31`, and `make_latex_tables.m:66` and `:126`.
- Modify: `docs/RUNNING_PIPELINE.md:80–150`.
- Test: `tests/PipelineStages.m:11–70`.

**Interfaces:**
- Consumes: `FilterByDualGraph(~curves)` (Task 5).
- Produces: the stages `FilterByDualGraphStar` (parallel, by level; input `curves_after_HHProposition1.dat`), `FilterByDualGraph` (parallel, by level; input `curves_after_UpdateCurves2.dat`) and `UpdateCurvesAfterDualGraph` (sequential closure).

- [ ] **Step 1: Write the failing test.** In `tests/PipelineStages.m`:
  1. In the `wc` list (lines 10–17), replace `"FilterByTraceStar", "HHProposition1", "FilterByTwistedTraceStar",` with `"FilterByTraceStar", "HHProposition1", "FilterByDualGraphStar", "FilterByTwistedTraceStar",`.
  2. In the same list, replace `"FilterByNonALInvolutionsStar", "UpdateByGenus", "UpdateCurves1",` with `"FilterByNonALInvolutionsStar", "UpdateByGenus", "UpdateCurves1", "FilterByALFixedPointsOnQuotient", "UpdateCurves2", "FilterByDualGraph", "UpdateCurvesAfterDualGraph", "Genus3CoversGenus2",`.
  3. In the `sh` list (lines 28–35), replace `"FilterByTraceStar", "HHProposition1", "FilterByTwistedTraceStar",` with `"FilterByTraceStar", "HHProposition1", "FilterByDualGraphStar", "FilterByTwistedTraceStar",`.
  4. In the same list, replace `"FilterByNonALInvolutionsStar", "UpdateCurves5", "FilterByAutomorphismGroup",` with `"FilterByNonALInvolutionsStar", "UpdateCurves2", "FilterByDualGraph", "UpdateCurvesAfterDualGraph", "Genus3CoversGenus2", "UpdateCurves5", "FilterByAutomorphismGroup",`.
  5. Change `assert #wstar eq 10;` (line 70) to `assert #wstar eq 11;`.

- [ ] **Step 2: Run it and watch it fail.** Run `magma -b filename:=tests/PipelineStages.m run_tests.m < /dev/null > /tmp/dg_task6.log 2>&1`. Expected: `PipelineStages.m: Fail!` with `Assertion failed`. The `workingcode.m positions` line shows a 0 for `FilterByDualGraphStar`.

- [ ] **Step 3: Edit `workingcode.m`** as listed under **Files**.

- [ ] **Step 4: Edit `run_pipeline.sh`.**
  1. In the header, after line 21 (the `writes input unchanged` continuation of `HHProposition1`), insert:
     `#   FilterByDualGraphStar                  [PARALLEL, by level]  (dual graph at p | D; D > 1, N squarefree)`
  2. In the header, after the `UpdateCurves2` line (35), insert:
     `#   FilterByDualGraph                      [PARALLEL, by level]  (dual graph at p | D)` and then
     `#   UpdateCurvesAfterDualGraph             [sequential]`
  3. After line 126 (`run_seq "HHProposition1" ...`), insert:
     ```bash
     # Dual graph of the special fibre at each p | D (Padurariu-Saia Thm 2.6; Baker-Norine/Chan):
     # no involution with tree quotient => not hyperelliptic.  After the check-only HHProposition1,
     # so VerifyHHTable2 still reads curves_after_FilterByTraceStar.dat.  Split by level.
     run_par "FilterByDualGraphStar"
     ```
  4. Replace lines 164 and following (`run_seq "Genus3CoversGenus2"  "${D}/curves_after_UpdateCurves2.dat" ...`) with:
     ```bash
     run_par "FilterByDualGraph"
     run_seq "UpdateCurvesAfterDualGraph" "${D}/curves_after_FilterByDualGraph.dat" \
                                                                                        "${D}/curves_after_UpdateCurvesAfterDualGraph.dat"
     run_seq "Genus3CoversGenus2"  "${D}/curves_after_UpdateCurvesAfterDualGraph.dat"   "${D}/curves_after_Genus3CoversGenus2.dat"
     ```

- [ ] **Step 5: Edit `run_parallel_filter.sh`.** In the `case "${STAGE}"`:
  1. Change the input of `FilterByTwistedTraceStar` (line 44) to `INPUT_DAT="${DATA_DIR}/curves_after_FilterByDualGraphStar.dat"`.
  2. Add these cases before `FilterByTwistedTraceStar)`:
     ```bash
         FilterByDualGraphStar)
             # Dual graph at p | D on the star curves; split by level (one fiber per (D, N, p)).
             INPUT_DAT="${DATA_DIR}/curves_after_HHProposition1.dat"
             ;;
         FilterByDualGraph)
             # Dual graph at p | D on all quotients, after UpdateCurves2; split by level.
             INPUT_DAT="${DATA_DIR}/curves_after_UpdateCurves2.dat"
             ;;
     ```
  3. In the `Supported:` echo lines, add `FilterByDualGraphStar, FilterByDualGraph,`.

- [ ] **Step 6: Edit `parallel_filter_worker.m`.**
  1. After line 58, add `star_of["FilterByDualGraph"]             := "FilterByDualGraphStar";`.
  2. In the level-grouped condition (lines 80–81), extend the set to
     `{"FilterByTwistedTrace", "FilterByTwistedWeilPolynomial", "FilterByTwistedTraceStar", "FilterByTwistedWeilPolynomialStar", "FilterByDualGraph", "FilterByDualGraphStar"}`.
  3. In the `case stage:` block (after line 112), add:
     ```magma
         when "FilterByDualGraph", "FilterByDualGraphStar":
             FilterByDualGraph(~subseq);
     ```
  4. In the header comment (lines 24–28), after "...the unit of work is a LEVEL", add the sentence: `FilterByDualGraph (and its Star version) is level-grouped too: the dual graph of X_0^D(N) at each p | D is built once per level.`

- [ ] **Step 7: Edit `run_sequential_stage.m`.** In the case at lines 62–66, add `"UpdateCurvesAfterDualGraph"` to the list of stages that call `UpdateCurves(~curves);`. The case then reads:
  ```magma
      when "UpdateCurves1", "UpdateCurves2", "UpdateCurves3", "UpdateCurves4",
           "UpdateCurves6", "UpdateCurves7", "UpdateCurves8",
           "UpdateCurvesAfterAutomorphismGroup", "UpdateCurvesAfterTwistedTrace",
           "UpdateCurvesAfterTwistedWeilPolynomial", "UpdateCurvesAfterDualGraph":
          UpdateCurves(~curves);
  ```

- [ ] **Step 8: Edit `Utils.m`.** In `CurveCostProxy`, after line 25 (`if g lt 3 then return R!0; end if;`), insert:
  ```magma
      if stage in {"FilterByDualGraph", "FilterByDualGraphStar"} then
          // Ideal classes of level N*p in the algebra of discriminant D/p (mass formula), times the
          // number of AL primes, summed over p | D: the class lookups dominate (spec section 8).
          // The worker takes the maximum over a level, whose fibers are shared.
          if not DualGraphApplicable(X`D, X`N) then return R!0; end if;
          return R!(#PrimeDivisors(DN) * &+[(&*[Integers() | q - 1 : q in PrimeDivisors(X`D div p)]) *
                                           (&*[Integers() | l + 1 : l in PrimeDivisors(X`N * p)])
                                           : p in PrimeDivisors(X`D)]);
      end if;
  ```

- [ ] **Step 9: Edit the stage lists.**
  1. `reconstruct_attribution.m`: after `"FilterByTraceStar",` (line 34) insert ` "FilterByDualGraphStar",`. After `"UpdateCurves2",` (line 48) insert ` "FilterByDualGraph",` and ` "UpdateCurvesAfterDualGraph",`.
  2. `analysis_stages.m`: after the HHProposition1 comment (line 20) insert `    "FilterByDualGraphStar",`. After `"UpdateCurves2",` (line 31) insert `    "FilterByDualGraph",` and `    "UpdateCurvesAfterDualGraph",`.
  3. `make_latex_tables.m`, in T1: after the `FilterByTraceStar` row (lines 65–66) insert:
     ```magma
      <"Dual graphs of the special fibres at $p \\mid D$ (Section~\\ref{sec:dualgraph})",
             "data/curves_after_FilterByDualGraphStar.dat">,
     ```
  4. `make_latex_tables.m`, in T2: after `<"Propagate closure and isomorphism", "data/curves_after_UpdateCurves2.dat">,` (line 126) insert:
     ```magma
      <"Dual graphs of the special fibres at $p \\mid D$ (Section~\\ref{sec:dualgraph})",
             "data/curves_after_FilterByDualGraph.dat">,
      <"Propagate closure and isomorphism", "data/curves_after_UpdateCurvesAfterDualGraph.dat">,
     ```
     `sec:dualgraph` is the label that the paper section for this test must carry. It states the unweighted variant only (user decision). Rows are dropped while their file does not exist (`make_latex_tables.m:80` and `:156`).

- [ ] **Step 10: Edit `docs/RUNNING_PIPELINE.md`.**
  1. In "Stage order", in the star block, after the `HHProposition1` line insert:
     `    FilterByDualGraphStar                     parallel, by level  NEW`
  2. In the all-quotients block, replace `    UpdateCurves2, Genus3CoversGenus2, UpdateCurves3` with:
     ```
         UpdateCurves2
         FilterByDualGraph                         parallel, by level  NEW
         UpdateCurvesAfterDualGraph                                    NEW
         Genus3CoversGenus2, UpdateCurves3
     ```
  3. After the paragraph "The star twisted stages decide some D = 1 curves...", add this paragraph:
     > `FilterByDualGraphStar` runs after the check-only `HHProposition1`, so `VerifyHHTable2` still reads `curves_after_FilterByTraceStar.dat`. The HH curves have D = 1 and are not applicable to it anyway. Both dual-graph stages apply only to D > 1 with N squarefree, and they skip everything else. They cache each fiber in `data/dualgraph/` (gitignored, keyed `dualgraph_v1_<D>_<N>_<p>`; override with `DUALGRAPH_CACHE_DIR`), so Phase B reuses the fibers built in Phase A. A curve whose quotient graph does not have the curve's genus stops the stage with `DUALGRAPH_GENUS`, like `BADDIM`.
  4. In "Splitting by level", change "The twisted stages (all four)" to "The twisted stages (all four) and the two dual-graph stages".
  5. In "Cost hot spots", add a bullet:
     > * **`FilterByDualGraphStar` / `FilterByDualGraph`**: one fiber per (D, N, p), cached on disk. The largest star level (39270, 1) takes 50 s at p = 2. The committed data gives about 2.7 CPU-h for Phase A (945 levels, 2300 fibers) and at most 1.2 CPU-h for Phase B before cache reuse. The per-curve work is under 0.01 s.

- [ ] **Step 11: Run the tests again.** Run the command from Step 2. Expected: `PipelineStages.m: Success!`. Then run `magma -b filename:=tests/DualGraph.m run_tests.m < /dev/null > /tmp/dg_task6b.log 2>&1`. Expected: `Success!`.
- [ ] **Step 12: Smoke test the worker.** Run the worker on one chunk of the committed star data:
  ```bash
  mkdir -p /tmp/dgsmoke && cp data/curves_after_UpdateByGenusStar.dat /tmp/dgsmoke/curves_after_HHProposition1.dat
  magma -b stage:=FilterByDualGraphStar chunk:=1 total_chunks:=512 input_dat:=/tmp/dgsmoke/curves_after_HHProposition1.dat \
        output_dat:=/tmp/dgsmoke/out.dat parallel_filter_worker.m < /dev/null > /tmp/dgsmoke/log 2>&1
  ```
  Expected in `/tmp/dgsmoke/log`: one `DualGraph summary:` line, and no `ERROR in worker`. `/tmp/dgsmoke/out.dat` exists. The chunk holds the most expensive levels, and was estimated at under 5 minutes.
- [ ] **Step 13: Propose a commit** with the message:
  `Pipeline: FilterByDualGraphStar after HHProposition1, FilterByDualGraph + UpdateCurvesAfterDualGraph after UpdateCurves2 (by level)`

## Task 7: SFI exclusion

**Files:**
- Modify: `ShimuraQuotients.m:1113`. After the paragraph ending `...were ruled at p = 3 through a source certified only by it.)`, add the comment.
- Test: `tests/DualGraph.m` (append).

**Interfaces:**
- Consumes: `NonHyperellipticAtPrimeCertificate(X::ShimuraQuot, p::RngIntElt) -> BoolElt, MonStgElt` (`ShimuraQuotients.m:1118`) and `FilterByDualGraph`.
- Produces: nothing new. This task adds a regression pin only.

- [ ] **Step 1: Write the test.** Append to `tests/DualGraph.m`:

<!-- file: tests/DualGraph.m append -->
```magma
procedure test_DualGraphNotSFICertificate()
    printf "Testing that DualGraph verdicts are not SFI certificates...";
    X := CreateShimuraQuot(158, 1, {1, 2}); X`g := GenusShimuraCurveQuotient(158, 1, {1, 2}); X`CurveID := 1;
    cs := [X];
    FilterByDualGraph(~cs : CacheDir := "none");
    assert cs[1]`TestInWhichProved eq "DualGraph at p = 2";
    // ruled by the graph at p = 2 | D, yet not certified at any good prime tried
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
```

- [ ] **Step 2: Run it.** Run `magma -b filename:=tests/DualGraph.m run_tests.m < /dev/null > /tmp/dg_task7.log 2>&1`. Expected: `...not SFI certificates...Done!` and `Success!`. This passes already, as a pin. To check that the pin has force, add `DualGraph` to a line inside the intrinsic body in a scratch copy, not in the tree, and confirm the text assert fails. Then discard the copy.
- [ ] **Step 3: Add the comment** to `ShimuraQuotients.m` after line 1113:
  ```magma
  // Deliberately NOT used: the dual-graph test (FilterByDualGraph, DualGraph.m).  It is a statement
  // about the special fibre at a prime p | D, where X has bad (Mumford) reduction, and says nothing
  // about X mod a good prime.  tests/DualGraph.m asserts this intrinsic never mentions it.
  ```
- [ ] **Step 4: Run the test again.** Expected: `Success!`. The comment sits above the intrinsic, so the text assert still passes.
- [ ] **Step 5: Propose a commit** with the message:
  `NonHyperellipticAtPrimeCertificate: document and pin that DualGraph verdicts are not admitted`

## Task 8: Stankewicz agreement, Brandt cross-check and the offline full-identity test

**Files:**
- Create: `tests/dualgraph_fiber_data.txt`
- Create: `tests/_offline/DualGraphFull.m`
- Test: `tests/DualGraph.m` (append).

**Interfaces:**
- Consumes: `DualGraphData`, `DualGraphQuotientFromData` and `ValidateDualGraphData`, plus Magma's `BrandtModule(D::RngIntElt, m::RngIntElt)` and `AtkinLehnerOperator(M::ModBrdt, q::RngIntElt)`.
- Produces: tests only.

- [ ] **Step 1: Generate the Stankewicz data.** Write this script to `/tmp/gen_fiber.m`:
  ```magma
  load "dual_graphs.m";
  out := [];
  for lev in [<6,35>, <10,33>, <14,5>] do
    D := lev[1]; N := lev[2];
    for p in PrimeDivisors(D) do
      ei, al := fiber(D, N, p);
      Append(~out, <D, N, p, [[x[1], x[2]] : x in ei], [<q, [a[2] : a in al[q]]> : q in Sort(Setseq(Keys(al)))]>);
    end for;
  end for;
  Write(OUT, "// Raw output of Stankewicz's fiber(D, N, p) (github.com/fsaia/GenusAtMost2, dual_graphs.m):\n// <D, N, p, [[origin, terminus] per edge], [<q, images of the edges under w_q>]>.\n// Generated by handoff_2026-09-29/graph-pilot/dual_graphs.m.\n" cat Sprint(out, "Magma") : Overwrite := true);
  quit;
  ```
  Run it from the pilot directory, because `load` is relative:
  ```bash
  cd /Users/sachihashimoto/Repos/ShimuraCurveALQuotients/handoff_2026-09-29/graph-pilot && \
    magma -b OUT:=$OLDPWD/tests/dualgraph_fiber_data.txt /tmp/gen_fiber.m < /dev/null > /tmp/gen_fiber.log 2>&1; cd -
  ```
  Expected: `tests/dualgraph_fiber_data.txt` is about 4.5 KB. The log shows harmless `AL_identifiers.m` / `quot_genus.m` load errors, which come from the pilot's `dual_graphs.m`. `fiber` does not need either file.

- [ ] **Step 2: Write the tests.** Append to `tests/DualGraph.m`:

<!-- file: tests/DualGraph.m append -->
```magma
dg_fiber_raw := eval Read("tests/dualgraph_fiber_data.txt");

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
        assert dg_same_graphs(DualGraphData(D, N, p : CacheDir := "none"), theirs, D, N, p);
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
```

- [ ] **Step 3: Create `tests/_offline/DualGraphFull.m`.**

<!-- file: tests/_offline/DualGraphFull.m create -->
```magma
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
```

- [ ] **Step 4: Run both.**
  * `magma -b filename:=tests/DualGraph.m run_tests.m < /dev/null > /tmp/dg_task8.log 2>&1`. Expected: all eight `Done!` lines and `Success!`, in about 16 s in total.
  * `magma -b filename:=tests/_offline/DualGraphFull.m run_tests.m < /dev/null > /tmp/dg_task8b.log 2>&1`. Expected: `DualGraphFull: 988 genus identities` and `Success!`.
- [ ] **Step 5: Run a negative control of the Stankewicz test.** In a scratch copy of `tests/dualgraph_fiber_data.txt` (not the tree), swap two AL images of one edge. Point a scratch copy of the test at that file, and confirm that it fails. Then discard both copies.
- [ ] **Step 6: Propose a commit** with the message:
  `tests: DualGraph agrees with Stankewicz's fiber and Magma's Brandt modules; offline 988 genus identities`

## Task 9: Handoff, plan pointer and the p = 2 count

**Files:**
- Modify: `HANDOFF.md:14`. Insert a new section above the first `## Handoff` heading.
- Modify: `PLAN.md:15`. Insert above `## ⇒ OPEN (2026-09-26)`.

**Interfaces:**
- Consumes: the stage names and the summary line of Task 5.
- Produces: documentation only.

- [ ] **Step 1: Add to `HANDOFF.md`** above line 14:
  ```markdown
  ## Handoff — 2026-09-29 — DUAL-GRAPH FILTER at p | D (FilterByDualGraphStar / FilterByDualGraph)

  * New package `DualGraph.m`. The dual graph of X_0^D(N) at p | D comes from left ideal classes
    (levels N and Np of the definite algebra of discriminant D/p). AL involutions are left
    multiplication by `PrimeIdeal`. The origin map is `lideal<O | Basis(I)>`, looked up by a
    canonical Gram key and confirmed by `IsIsomorphic`. It is unweighted only.
  * Validation: 988/988 genus identities (`tests/_offline/DualGraphFull.m`; 822 in CI), K4 for
    X_0^6(55)/<w_3, w_110> at p = 3, ten hyperelliptic controls, three negative controls,
    agreement with Stankewicz's `fiber` (300 quotient graphs), and Brandt fixed-point counts.
  * `BrandtModule(D', M)` exposes no ideals and no degeneracy maps, and its basis cannot be aligned
    with `LeftIdealClasses`. It is also slow: `BrandtModule(3,1610)` takes 232 s, against 10 s for
    our whole fiber. It is used as a cross-check only.
  * `LeftIdealClasses` labels differ between calls: compare builds up to isomorphism.
  * Verdicts are tagged `DualGraph at p = <p>`, with odd primes tried first. Count the verdicts that
    rely on p = 2 alone with:
    `c := eval Read("data/par/curves_after_UpdateCurves8.dat"); #[X : X in c | assigned X`TestInWhichProved and X`TestInWhichProved eq "DualGraph at p = 2"];`
    The pilot measured 35 of 165.
  * NOT an SFI certificate (pinned in `tests/DualGraph.m`).
  * Cache: `data/dualgraph/dualgraph_v1_<D>_<N>_<p>`, gitignored and regenerated on lava. Bump
    `DG_FORMAT` if the construction changes.
  * Next: full rerun on lava. Paper section `sec:dualgraph` (unweighted statement). Paper totals
    in the later pass.
  ```
- [ ] **Step 2: Add to `PLAN.md`** above line 16 (`## ⇒ OPEN (2026-09-26): ...`):
  ```markdown
  ## ⇒ OPEN (2026-09-29): dual-graph filter awaiting the lava rerun

  Implemented on branch `dual-graph-filter` (from `integration`). Plan:
  `docs/superpowers/plans/2026-09-29-dual-graph-filter.md`; spec:
  `docs/superpowers/specs/2026-09-29-dual-graph-filter-design.md`. To do: rerun on lava into an
  empty data dir. Report the `DualGraph summary:` lines of both stages and the p = 2-only count
  (command in HANDOFF.md, 2026-09-29). Then write the paper section `sec:dualgraph`.
  ```
- [ ] **Step 3: Check the text.** Run `grep -n "DualGraph" HANDOFF.md PLAN.md`. Expected: the lines above, and nothing else new.
- [ ] **Step 4: Propose a commit** with the message:
  `HANDOFF/PLAN: dual-graph filter state and the p = 2 count`
- [ ] **Step 5: Stop.** Ask the user whether to open a PR from `dual-graph-filter` into `integration`. Eran merges. Do not push without a yes.
