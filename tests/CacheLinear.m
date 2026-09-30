// The store-backed caches (Caching.m SetCache, ClassNumberData.m clRemember and
// ClassNumberBatchLU, GeneralizedComplicatedFixedPoints.m, CollectDisc) used to mutate an
// associative array the store still referenced.  Magma copies a value on write when more than
// one reference exists, so every insert copied the whole array: quadratic in cache size.  The
// fix takes the array out of the store before the insert.
//
// Sources for the numbers here: a Mac (Apple silicon, Magma 2.29), 2026-09-30, scripts
// vvdata/weyl-campaign/weil-cost-2026-09-30/cache_bench{,2}.m on the campaign branch.  Old
// SetCache: 20k inserts 4.3 s, 40k 17.2 s, 80k 68.2 s (a plain associative array: 0.01 s).
// Fixed: 80k in 0.15 s, 1.28M in 2.4 s.
//
// (1) uses a PRIVATE store, never the shared class_nos/CL_ORDER caches, so a failure here
//     cannot leave wrong values for later tests in the same run.
// (1) 80k inserts finish in under 10 s: 60x the fixed cost measured above, 6x below the old
//     one, so it separates the two behaviours on any machine this suite runs on.
// (2) cache semantics on that private store.
// (3) ClassNumberLU's own cache: the values agree with ClassNumber, and after one pass every
//     discriminant is a KEY of its cache (checked directly, not by timing a second pass, which
//     would be too fast to distinguish a hit from a recomputation).

_ := ClassNumberLU(-4);
import "Caching.m" : SetCache, GetCache;
import "ClassNumberData.m" : clGetAssoc, CL_ORDER;

// (1)
st := NewStore();
t0 := Cputime();
for i in [1..80000] do SetCache(-4*i, i, st); end for;
t := Cputime(t0);
printf "  80000 SetCache inserts into a private store: %o s\n", RealField(3)!t;
assert t lt 10;

// (2)
for i in [1, 2, 40000, 80000] do
    b, v := GetCache(-4*i, st);
    assert b and v eq i;
end for;
b, _ := GetCache(-4*80001, st);
assert not b;
SetCache(-4, 999, st);
b, v := GetCache(-4, st);
assert b and v eq 999;
CacheClear(st);
b, _ := GetCache(-4, st);
assert not b;

// (3) -- on a cleared cache, so an earlier test in the same run cannot have filled it.
CacheClear(CL_ORDER);
ds := [-4*k : k in [1..300]] cat [-4*k - 3 : k in [1..300]];   // -4k (0 mod 4) and -(4k+3) (1 mod 4)
assert &and[d mod 4 in [0, 1] : d in ds];
h := [ClassNumberLU(d) : d in ds];
assert h eq [ClassNumber(d) : d in ds];
order := clGetAssoc(CL_ORDER);
assert &and[IsDefined(order, d) : d in ds];
assert [order[d] : d in ds] eq h;
printf "  %o ClassNumberLU values agree with ClassNumber and are all cached\n", #ds;
