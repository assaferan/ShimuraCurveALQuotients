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

Those **23 will each pay a full Normaliz solve** rather than a cache hit — and, per
[[polymake-breakage-normaliz-fix]], an uncached level *silently* costs that solve instead of
erroring. Worth pre-solving them, or at least expecting them to dominate the wall-clock, before
wave 4b is launched at the original 2400 s cap.
