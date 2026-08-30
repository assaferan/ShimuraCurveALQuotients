# `MISSING_TARGETS.txt` — the broad target list, and what it says about the cache

Rescued 2026-08-30 from the `-diagnostic` worktree, where it had sat untracked since
2026-07-30. It was the **only copy** and is not superseded by the wave-4 lists.

Format is `base genus M`, one per line, where **`M = D·N`**:

    6_25 7 150        # D=6, N=25, genus 7, M = 150

## It is broader than the wave-4 lists, not a subset

    MISSING_TARGETS.txt      351 bases
    wave4_* (all lists)      140 bases
    overlap                   46 bases

So this is a wider sweep than anything the wave-4 triage covered — closer to the 377/73/304
counts in the post-m0 backlog plan than to the 140 the campaign has actually been driving.
Treat it as the superset to triage against, not as a stale duplicate.

## Cache implication for wave 4b

The `M` column is exactly the quantity the Normaliz solution cache is indexed by, and the
committed frontier on `main` is **M ≤ 2260** (`f90c441`, plus the two M = 532 files added in
`9cc771e`). Against that frontier:

    within the frontier (M <= 2260)    328 of 351
    beyond it (M > 2260)                23 of 351     max M = 4830

Those 23 pay a full Normaliz solve rather than a cache hit — and, per
[[polymake-breakage-normaliz-fix]], an uncached level *silently* costs that solve instead of
erroring. Worth pre-solving before wave 4b is launched at the original 2400 s cap.

### It is 11 pre-solves, not 23

The cache is keyed by `M`, which is the first component of `polymake_solution_<M>_<a>_<b>`, so
one solve serves every base at that level. The 23 bases collapse to **11 distinct M**, and the
sharing is heavy:

    2262  2262_1
    2310  210_11  330_7  462_5  6_385  770_3      <- 5 bases, one level
    2370  2370_1
    2478  2478_1
    2490  2490_1
    2670  2670_1
    2730  210_13  390_7  546_5  910_3             <- 4 bases
    3570  210_17  510_7  714_5                    <- 3
    3990  210_19  570_7  798_5                    <- 3
    4290  330_13  390_11                          <- 2
    4830  210_23

Cached levels carry 2–3 solution files each (M = 2260 → 2, M = 2076 → 3), so the real cost is
roughly **22–33 solves covering all 23 bases** — a bounded, batchable job rather than a per-base
tax, and one that can be run ahead of wave 4b instead of inside it.

### These are the high-genus bases, so expect CM starvation too

Genus distribution of the 23: g = 9 (×11), 11 (×4), 13 (×3), 15 (×4), 17 (×1). Per
[[cm-supply-second-triage-axis]] the CM-point demand is `max(2g+5)` over covers, not the MaxNum
default of 7 — at g = 17 that is 39. So even with the cache warm, this cohort is the one most
likely to die in "Could not find enough points" rather than on polytope time. Run `cmsupply.m`
over these 23 before committing wall-clock to them; `Targets` is the rescue lever.
