# Rescued probes — two patches with caveats

Both files below were **single-copy in throwaway worktrees** and are preserved here for the
reasoning, not as things to apply blindly. Read the caveats before using either.

## `bfprof-instrumentation.patch` + `bfprof.m` + `dsmall.m`

Rescued 2026-08-30 from the `-spanprobe` worktree, which `HANDOFF.md` marks THROWAWAY.

`bfprof.m` is the per-stage profiling driver that produced the odd-D profile — the one that
attributed 93% of an `m_idx` pass at `65_2` to the single unguarded recombination at
`BorcherdsForms.m:817`. It wraps `BorcherdsForms` in a `try`/`catch` so a runtime error still
prints the `SUMMARY` line; `magma -b` otherwise hangs on an uncaught error. It is the template
named by the handoff for the **invariant-hoist** measurement (the 336× redundancy in the 0-side
block), so it needed to survive the worktree.

    BFPROF=1 magma DD:=65 NN:=2 bfprof.m

`dsmall.m` prints `#disc_grp` for eight small bases — the sanity table for `weildim2.m`'s O(d)
trace formulas.

**CAVEAT — the patch does NOT apply to current `main` and is a record only.** It was taken
against the `-spanprobe` tree (based on `a924e1a`) and diffed against `origin/main`, so its 18
hunks mix two unrelated things:

1. the **`BFPROF`/`BFINV` timers**, which is the part worth reusing; and
2. an inline **prototype of the `WeaklyHolomorphicBasis` speedup** (`_p := 1073741789`,
   `_sel`, `_Csel`, …) that was **superseded** by the version merged to main as `04f1d7b`.

Applying it whole would duplicate the speedup. Lift the timers by hand.

## `valuesatcmpoints-characterization.patch`

Rescued 2026-08-30 from the `-diagnostic` worktree (branch `valuesatcmpoints-diagnostic`),
where it was uncommitted.

Turns the first-failure `error` in `ValuesAtCMPoints` into a **characterization pass**: collect
*every* non-rational cell rather than dying on the first, then correlate against
`gcd(disc, N)` and `KroneckerSymbol(d, p)` at the primes of `N`. The intent was to test whether
non-rationality tracks ramification at the level prime — see [[embedding-selection-root-cause]].

**CAVEAT — this probe cannot have run successfully as written.** It patches the two-argument
intrinsic

    ValuesAtCMPoints(abs_schofer_tab::SchoferTable, all_cm_pts::SeqEnum : Exclude := {})

at `SchoferFormula.m:1498`, which has **no `Xstar` in scope** — that is the *four*-argument
intrinsic at `:1655`. The added lines reference `Xstar`N` and `Xstar`D`, so the patch fails on
an undefined identifier the moment a bad cell is found. **Any conclusion attributed to this
probe should be treated as unevidenced.** To revive it, thread `D` and `N` through the
signature, or move the characterization up to the four-argument caller.
