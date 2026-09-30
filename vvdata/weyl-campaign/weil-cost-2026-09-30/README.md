# Weil-stage cost measurement, 2026-09-30

Where `FilterByWeilPolynomial` spends its time, replayed per prime exactly as the stage does
(`weil_timing.m`: `WeilPolynomial(X, p)` over `PrimesUpTo(bound)`, bound = min(table, budget,
ceiling); each prime's CPU/wall time written with `Write()` so it survives stdout buffering).

    magma -b id:=<CurveID> tag:=<label> scratch:=<logdir> weil_timing.m     # from a code tree

## Why these numbers are trustworthy, and what they do NOT say

* Same Weil polynomials on every tree (`main`, #56, #58) -- the runs cross-check each other.
* ⚠ Taken on a Mac with `CLASS_GROUPS_FAST_DIR` UNSET: every new discriminant was a direct
  `ClassNumber` (21770 of them = 20 s in one n = 13^6 call, `prof_weil_leaf.log`).  Pipeline runs
  happen on lovelace/lava with the tables mounted.  These numbers support the ORDERING of curves
  (which shape is heavy) and the RATIOS between code trees on the same curve; they are not the
  stage's absolute times.  Re-measure on lovelace before quoting an hour figure.

## Results (CPU s, `main` @ 95cf87b unless labelled)

    curve  shape                              main               #56 (memo)          #58 (cache fix)
    13029  X_0^210(73)/W32  g=4 Qmax=15330    559 (7 primes)     273                 --
    1071   X_0(240)/W4      g=6 Qmax=80       2959 thru p=17     704 thru 17         1171 thru 17; 2153 thru 19
                                              (~9 h projected)   (p=19: 1479)        (p=7/11/13/17/19: 41/166/288/677/982)
    7296   X_0^21(20)/W4    g=7 Qmax=20       --                 p=11: 201, p=13: 1137

    per-prime, main, curve 1071:  p=7 43, p=11 224, p=13 560, p=17 2132
    per-prime, memo, curve 1071:  p=7 3.8, p=11 35.6, p=13 119, p=17 546, p=19 1479

Cost is set by n = p^g (n=9e5: 206 s, 4.8e6: 560 s, 2.4e7: 2132 s on main, across curves),
nearly independent of Qmax / W / level; exponent in n between 0.8 and 1.5 depending on the curve.

## What came out of it

* PR #58: the store-backed caches copied the whole associative array on every insert
  (`cache_bench.log`: 80k SetCache inserts 68 s; `cache_bench2.log`: 1.28M in 2.4 s once fixed).
* PR #59: `CurveCostProxy`'s Weil branch summed 4*Qmax*p^g (the discriminant DEPTH) and ranked
  the 9-min curve above the 9-h one, 1071 at rank 138/886 (`proxy_check.log`); sum p^g puts
  7296/1071 at ranks 1/2 (`proxy_alt.log` used p^(0.8 g); the PR uses p^g).
* The `docs/RUNNING_PIPELINE.md` "about 5 h on a single curve" is the right order for the
  genus-6 small-W shape, not the big-level shape it named.
