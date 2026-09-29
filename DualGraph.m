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
