# `IntegralSolution`: the NONINTEGRAL class, measured

Result of making the 2026-08-29 `intsol` finding usable. Code landed on branch `intsol-optin`
(`6a6267c`, off `main`); the reverted prototype it descends from is
`../intsol-acceptance-criterion.patch`.

## What was actually wrong before

The finding was never in doubt — what blocked it was that it shipped **unconditionally**. It
replaces the solve's solution, so `fs[-1]` becomes a different form and the reference
comparisons in `tests/SchoferIsometry.m` (Guo–Yang Table 45) and `tests/VectorValuedForm.m`
(the 15_2 multiplier) fail, on bases that are already fine. The prototype's own comment says why
it was unconditional: **`assigned` cannot see a top-level variable from inside a package
intrinsic**, so a flag would never have fired. A named **parameter** threads correctly. That is
the entire fix.

Default `false` ⇒ inert, bit-for-bit unchanged. All three tests that call `BorcherdsForms` pass
(`SchoferIsometry` 25.5 s, `M0MultiplierExact` 18.6 s, `VectorValuedForm` 97.3 s). Threaded
through `EquationsOfCovers` and both `AllEquationsAboveCovers` forms so model generation can opt
in (`INTSOL=1` in `genmodels.m`).

## The measurement (`intsol.m`, both passes per base)

    RESCUED (integral, den 1)   33_2  34_11  46_3  57_2  58_5  74_3  74_5
    better, not cleared         14_37 (6->2)   39_2 (65536->32)   82_3 (3->2)
    unchanged                   22_17  46_5  51_2  62_3  62_7  86_3
    WORSE                       69_2 (64 -> 1048576)
    rejected every triple       26_7
    no verdict                  111_2 (capped at 2h32m)   87_2 (still running)

**≈39% of the class is a choice artifact** — the non-integrality was which point `Solution`
returned from `sol + Kernel`, not a property of the divisor.

## Two independent confirmations of earlier work

* **`74_3` is RESCUED.** [[nonintegrality-origin-echelonform]] used that very base to conclude
  the defect is in the SOLUTION the divisor solve picks, not the basis. Fixing the solution
  rescues it — exactly as predicted.
* **`39_2` improves but does NOT clear** (65536 → 32). [[nonintegral-forms-root-cause]] records
  `39_2` as a **basis-level** failure `intsol` cannot touch. Also exactly right.

Both were written from different evidence, months apart, and both hold.

## THE CAVEAT: it is not monotone — keep it opt-in

**`69_2` gets four orders of magnitude WORSE** (denominator 64 → 1048576). An integral *solution
vector* against a *non-integral basis* can produce worse principal parts than the arbitrary
rational solution did — the two objects are not the same, and `EchelonForm` is what makes the
pool non-integral ([[nonintegrality-origin-echelonform]], denominators ~3e28; `111_2`'s default
pass reports ~3.2e30, which is that pathology rather than ordinary arithmetic).

So **do not promote this to an automatic fallback.** Per-base opt-in is the right granularity.

`26_7` shows the other edge: rejecting non-integral triples can leave *nothing*, moving a base
from NONINTEGRAL to form-failure rather than to a model. Honest, but not obviously an
improvement.

## Necessary, not sufficient

Integral principal parts do not make a model. The seven still have to clear CM supply and the
`y2`-scale assembly — and `58_5` sits right beside the documented CM-starved `58_3`
([[cm-supply-second-triage-axis]]). Anything that does produce a model must go through
`ModelVerification` (genus, Weil-polynomial divisibility, trace-formula point counts), none of
which touches the Borcherds/Schofer path that produced it.

## Method note

The concurrency limit in the first sweep runner (`zsh`, `jobs -r | wc -l`) **silently did not
work** in a non-interactive script — all 16 jobs launched at once (load 33 on 14 cores). Memory
was fine and the verdicts are unaffected, because this sweep measures integrality and
denominators rather than durations — but **the timing column in these logs is confounded and
must not be quoted.** Use `xargs -P N` for a limit that actually holds.

## Files

    intsol.m           the two-pass driver (default vs IntegralSolution)
    sweep_<base>.log   all runs
