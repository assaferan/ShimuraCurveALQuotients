# Route A is NOT the cheap predictor for the 81 unclassified bases — measured

The 122-base sweep left **81 bases with no verdict** because `ppint` must construct Borcherds
forms before it can speak, and on large bases that exceeds an hour. The obvious escape was
route A (`deficit.m`): the Borcherds-obstruction deficit is
`Ncols(mat) - Rank(ech_basis * mat)`, which depends only on the weakly holomorphic basis and
`coeffs_to_divisor_matrix` — **no CM points, no divisor-triple loop**.

[[deficit-predictor-two-routes]] recorded route A's whole cost as *"one `EchelonForm` of a
12784x258 matrix = 755 s of 810 s"*, and the 66x `WeaklyHolomorphicBasis` speedup (`04f1d7b`)
took exactly that echelon from **755 s to 1.5 s**. So the predictor looked like it should have
become cheap for free.

**It did not.** Measured on three bases that produced no `ppint` verdict in 3600 s:

    10_67    WeaklyHolomorphicBasis  1242 s   (n = 565, n0 = 101, #eta = 467)
    10_71    WeaklyHolomorphicBasis  1430 s   (n = 599, n0 = 121, #eta = 495)
    10_101   capped at 1800 s -- did not finish WeaklyHolomorphicBasis at all

`WeaklyHolomorphicBasis` still **dominates** at 20-24 minutes. The speedup fixed the echelon,
but on this cohort the cost simply lives elsewhere: these bases carry `#eta ~ 500` and pole
order `n ~ 600`, far beyond `38_5` (where the 66x was measured). The script's own header
anticipated exactly this — *"if `WeaklyHolomorphicBasis` dominates, the predictor costs what the
pipeline costs"* — and it does.

**⇒ Route A does not rescue the 81.** A cheap predictor for this cohort must avoid building the
basis at all, which leaves route B (the Weil-representation dimension formula). Per
[[deficit-predictor-two-routes]] route B's `k = 3/2` phase is wrong by nearly exactly `-d/6`,
with an explicit warning **not to tune that constant to fit**.

## SEPARATELY: `deficit.m` is BROKEN on the campaign branch

    User error: Identifier 'ProbeWHBasis' has not been declared or assigned   (deficit.m:50)

It calls `ProbeWHBasis` and `ProbeDivisorMatrix`, which existed **only in the `-spanprobe`
worktree's patched `BorcherdsForms.m`**. That worktree was retired on 2026-08-30 (its content
preserved as `bfprof-instrumentation.patch`), so `deficit.m` has been non-functional since.

This matters twice over:

* the three validated results in [[deficit-predictor-two-routes]] (`38_5` -> 1, `38_7` -> 0,
  `34_3` -> 0) were obtained **with that patch applied**, which the note does not say;
* anyone running `deficit.m` today pays the **full `WeaklyHolomorphicBasis` cost first** — 20+
  minutes on a large base — and only then hits the undefined-identifier error.

To revive it, either apply `bfprof-instrumentation.patch` or replace the two `Probe*` helpers
with their production equivalents (`basis_of_weakly_holomorphic_forms`, `coeffs_to_divisor_matrix`),
which is the smaller change and does not depend on a deleted worktree.
