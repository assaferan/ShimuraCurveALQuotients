# Gate 4 (class-constancy) at 34_11: a GENUINE violation, and it tracks the CLASS

`X_0(34,11)` clears gates 1-3 and then fails

    M0MultiplierExact: class-constancy violated (class 1, dev 0.0185658, scale 0.0430097)

**43% of scale.** The check's own comment records the empirical roundoff floor as `1.1e-22` at
the deepest known base (`39_2`, M = 156), whose independent dump passes every class-constancy
check (42/42) — so constancy itself is exact — and states outright that *a genuine violation is
O(scale)*. This is twenty orders past that floor. **Loosening the guard would manufacture a
multiplier wrong by 43%.**

## What the check compares

`contribs[k][j] = rvtab[wi][i] * a0w[wi]` — the contribution of coset `wi` to isotropic component
`j`. Theory says this is constant across all cosets of a cusp class (the result
[[cusp-class-assembly-closed]] established at 570/570). The code samples several cosets per class
(`picks`) and compares each against the first.

## The measurement

Reporting every comparison AND still erroring, so the printed data covers exactly the
comparisons made before the failure fires:

    class 1, picks 1/2/3 -- three sampled cosets of ONE class
    pick 1 (wi 2):  |contrib| 0.043009714   |rv| 0.0018906598   |a0w| 22.748520
    max dev across picks: 0.0185658   =  43% of scale

**Three cosets the code believes share a cusp class produce materially different
contributions.** The deviation tracks the CLASS, not one rogue coset — which points at the
**cusp-class partition (`classes` / `canon`) being wrong for this base**, rather than at
`rvtab` or `a0w` being miscomputed at a particular coset.

Not yet proven: confirming that neither factor is individually the outlier needs `|rv|` and
`|a0w|` compared across picks 2 and 3. That is a small extension of the same probe.

## A METHOD WARNING, recorded because it produced a false negative

A first pass ran the probe in **report-and-continue** mode and showed **max dev ~0 everywhere** —
flatly contradicting the observed failure, and I briefly reported gate 4 as possibly unreal on
that basis. Two separate causes:

* the analysis read the wrong awk field (`$12` is the literal string `dev`, not its value), so
  everything parsed as zero;
* and reporting *instead of* erroring changes which comparisons are reached.

**Run the probe in report-AND-error mode**, so the printed set is exactly the set the production
check evaluates. The two views must be made to coincide before either is trusted.

## The probes themselves

`gate3b-gate4-probes.patch` (against `intsol-optin` @ `969fa85`) carries both diagnostics:

* **`GATE3B=1`** — report every failing slash-constant comparison and continue, instead of
  erroring. Used to ask whether `10_61`/`14_43` fail gate 3 for a reason a *third* tolerance
  change could fix. Production path is untouched when the variable is unset.
* **`GATE4=1|other`** — report class-constancy comparisons. `GATE4=1` reports and continues;
  any other value reports AND still errors, which is the mode that made the two views coincide.

Both are THROWAWAY instrumentation and must not be merged. They exist because
report-instead-of-error cracked gate 3, exposed the lost logs in the local sweep, and settled
gate 4 — but the same technique produced a false negative when used carelessly (see the method
warning above).
