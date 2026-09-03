# Why `cusp7_15_2.out` has no `g=2,6,10,30`: not a dump bug, a monomial-pool gap

`PLAN.md` (tier1-models, MAIN LINE) diagnosed the missing classes as `cusp7.m`'s
"first-encountered-per-class" `ppdumped` logic: `dopp` was set once per class from the
*first* coset word seen, so if that word had `lead > 0` (no pole/constant term) for some
monomial, later words in the same class could never backfill it.

That bug is real and is now fixed (`ppdumped` keyed on `<cls0, ri>` instead of `cls0`
alone) — but **it is not why classes 2, 6, 10, 30 are empty.** Re-running the fixed
`cusp7.m` on `15_2` produces a byte-identical `PP` section (3557 lines, same content;
the only diff is floating-point summation-order noise at ~1e-71 in the `SELFC`
self-check, six orders below anything that matters).

**Directly verified cause:** for *every* coset word whose `cls0 = GCD(g[2][1] mod M, M)`
is 2, 6, 10, or 30, *every one of the 158 monomials* (9-form panel ∪ 69 hand-picked
`extraseqs`) has `lead > 0` — i.e. strictly holomorphic there, no pole, no constant term
to record. Measured (`M = 60`, 15_2):

    cls0    words         W    minlead   #monomials with lead<=0
    2       62,63,...,76  30   36        0
    6       124..128      10   12        0
    10      129..131       6   36        0
    30      144            2   12        0

All 12 divisors of `M = 60` genuinely occur as `cls0` among the 144 coset words (checked
directly) — the classes exist, the current monomial pool just never touches them.

**Consequence for the `rem:gauge` resolve:** the "concrete next step" in `PLAN.md` —
"re-run a `cusp7.m`-style pass that guarantees every intermediate class gets dumped" — is
not enough on its own; the dedup fix here is necessary but not sufficient. What's
actually needed is new eta-monomials (valid weight-appropriate exponent vectors, same
12-component shape as the current pool) that carry a pole or constant term at cusps
2, 6, 10 and 30 specifically. `extraseqs`'s 69 vectors were hand-picked in an earlier,
unrecorded session with no surviving generator script — reconstructing or extending that
pool is itself the open task, not a data-completeness bug. Scope it as a real search
(likely: query a broader weakly-holomorphic-basis pool at these four cusps, or derive
which exponent-vector constraints put a pole exactly at a given cusp class from the
`triang`/`leads` formula in `cusp7.m` directly) before attempting it.

The `ppdumped` fix is kept anyway: it is strictly more correct in general (protects any
future base/monomial-pool combination where a first-encountered word legitimately misses
some monomials that a later word in the same class would catch), even though it happens
to be a no-op on this specific 15_2 run.
