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
are involutions of the right sets, w_p commutes with every other w_q on edges, the origin map is
equivariant for every w_q with q | DN/p, and h' - 2h + 1 is the genus of X_0^D(N).  Run on every
build and on every disk-cache read.}
    require #data eq 5 : "data must be <h, h', org, alE, alV>";
    h := data[1]; hh := data[2]; org := data[3];
    assert #org eq hh and &and[v in [1..h] : v in org];
    assert {x[1] : x in data[4]} eq Set(PrimeDivisors(D*N));
    assert {x[1] : x in data[5]} eq Set(PrimeDivisors((D div p)*N));
    wp := dg_img(data[4], p);
    for x in data[4] do
        img := x[2];
        assert Sort(img) eq [1..hh] and &and[img[img[e]] eq e : e in [1..hh]];
        // w_p commutes with every w_q on edges: the vertex action on the - copy uses alV[q]
        assert &and[wp[img[e]] eq img[wp[e]] : e in [1..hh]];
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
of Hall divisors of DN, a group, so containing 1), with half-edges removed and leaves kept.  Returns
the number of vertex orbits that carry an edge and the edges as pairs <u, v>, u <= v.}
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
