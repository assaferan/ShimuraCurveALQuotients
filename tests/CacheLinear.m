// The store-backed caches (Caching.m SetCache, ClassNumberData.m clRemember / clFundClassNo,
// GeneralizedComplicatedFixedPoints.m) used to mutate an associative array the store still
// referenced.  Magma copies a value on write when more than one reference exists, so every
// insert copied the whole array: quadratic in cache size.  Measured on main at 95cf87b,
// SetCache alone: 20k inserts 4.3 s, 40k 17.2 s, 80k 68.2 s.  The class-number cache reaches
// hundreds of thousands of entries over a Weil-stage run, so this, not the trace formula,
// was where that stage's hours went.  Fix: drop the store's reference before mutating.
//
// (1) SetCache stays linear: 80k inserts well under the old 68 s.  The threshold is 60x
//     above the fixed cost measured here (0.15 s) and 6x below the old cost, so it separates
//     the two behaviours on any machine this suite runs on.
// (2) The cache is still correct: every key readable, overwrite works, misses miss, and a
//     CacheClear empties it.
// (3) ClassNumberLU's own cache: repeated lookups are hits (no recomputation) and values agree
//     with ClassNumber.

_ := ClassNumberLU(-4);
import "Caching.m" : SetCache, GetCache, class_nos;

// (1)
CacheClear(class_nos);
t0 := Cputime();
for i in [1..80000] do SetCache(-4*i, i, class_nos); end for;
t := Cputime(t0);
printf "  80000 SetCache inserts: %o s\n", RealField(3)!t;
assert t lt 10;

// (2)
for i in [1, 2, 40000, 80000] do
    b, v := GetCache(-4*i, class_nos);
    assert b and v eq i;
end for;
b, _ := GetCache(-4*80001, class_nos);
assert not b;
SetCache(-4, 999, class_nos);
b, v := GetCache(-4, class_nos);
assert b and v eq 999;
CacheClear(class_nos);
b, _ := GetCache(-4, class_nos);
assert not b;

// (3)
ds := [-4*k : k in [1..300]] cat [-4*k - 3 : k in [1..300]];   // -4k (0 mod 4) and -(4k+3) (1 mod 4)
assert &and[d mod 4 in [0, 1] : d in ds];
h1 := [ClassNumberLU(d) : d in ds];
t0 := Cputime();
h2 := [ClassNumberLU(d) : d in ds];
t := Cputime(t0);
assert h1 eq h2;
assert h1 eq [ClassNumber(d) : d in ds];
printf "  %o ClassNumberLU values agree with ClassNumber; repeat pass %o s\n", #ds, RealField(3)!t;
assert t lt 1;
