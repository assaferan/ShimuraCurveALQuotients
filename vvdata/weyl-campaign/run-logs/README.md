# Long-run logs kept as provenance

Runs too long to reproduce casually, kept so a future session can check what a recorded number
actually came from rather than re-spending the wall-clock.

## `x0_111_1-lava-2026-09-12.log` — 58595 s (16.3 h) on lava, PASSES

The second of the two Guo-Yang bases that had no re-derivation test (`93_1` was the first). It
anchors on the FULL genus-7 hyperelliptic curve, which is a stronger anchor than `93_1`'s
quotient-only one:

    X0^111(1): 1 curve comparison(s), 0 involution comparison(s), 1/1 expected covers matched;
               4 committed model cover(s) re-derived (0 CRV skipped)

⚠ **It survived the pool wall.** The handoff flagged its pool reaching 1678 vectors at `m_idx=3 of
7` against a recorded ~2000-vector / ~11 GB death point. It did not die -- so that wall is a
recorded observation, not a hard limit at this size.

⚠ **The sign fallback fired and recovered**: `sign of -88 unresolved (2 candidates); trying spare
-168`. Worth knowing if this base is ever re-run and the spare discriminant changes -- the run is
not sign-deterministic in the way a base with no ambiguity would be.
