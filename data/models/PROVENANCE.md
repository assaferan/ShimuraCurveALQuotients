# How to regenerate every committed model

Single source of truth for **which flags each model file needs**. Written 2026-09-05, after a
session that found five files unregenerable under default settings — and, worse, no record
anywhere of what they *did* need.

**The invariant this file exists to protect:** every file in `data/models/` must be regenerable by
a command written down here. If you produce a model with a non-default flag, add its row before
committing the model.

## The default recipe

    export NORMALIZ_BIN=~/Documents/GitHub/normaliz-3.11.1/normaliz
    magma -b D_s:=<D> N_s:=<N> OUTDIR:=<dir> genmodels.m < /dev/null

`genmodels.m` lives on the `m0-theta-campaign` branch at `vvdata/weyl-campaign/genmodels.m`
(triage tooling is not on `main` — see `CLAUDE.md`). **Always redirect stdin**, or a runtime error
drops Magma into its interpreter and blocks forever.

## Files needing a non-default flag

| base | flag | why | validated by |
|---|---|---|---|
| `39_2` | `CMNONCOPRIME=1` | default sees 3 CM points against demand 19 and dies with "Could not find enough points"; with the coprime-to-level filter off it sees 24 and builds 15 keys, 0 empty | `tests/GuoYangEquations.m` — `W={1}` genus-7 curve `IsIsomorphic` to Guo-Yang's published equation |
| `14_3` | `CMNONCOPRIME=1` | default leaves the covers under-determined and `W={1}` EMPTY (6 keys, 3 populated); with the filter off, 16 keys, 0 empty | `tests/_offline/GuoYangCurve_14_3.m` — `IsIsomorphic` to Guo-Yang, exact, 6817 s |
| `22_5` | `Y2TWIST=1` | the unpinned-y2-scale guard (`1768517`) POSTDATES the file, so default regeneration withholds `[1,2,5,10]` | restores 3/3 covers, coefficient-for-coefficient |
| `15_2` | `Y2TWIST=1` | same guard | restores 12/12 keys (one cover up to isomorphism) |
| `22_3` | `Y2TWIST=1` | same guard | restores 13/14; `[1,3,22,66]` and one `[1,66]` entry are **still lost** — both genus-0 conics, which `select_y2_twist` skips by construction |
| `26_3` | `CMNONCOPRIME=1` | the default filter admits only 3 of Guo-Yang's 14 discriminants (only `-8, -11, -20` are coprime to `N=3`) against demand 9, and the base dies at "Could not find enough points"; with the filter off the pool is 21 | `tests/GuoYangEquations.m` — the full `V_4` diagram matches Guo-Yang: `[1,78]` and `[1,3]` `IsIsomorphic` to their genus-2 and genus-3 quotients, and `[1,26]` entry 2 IS their conic `-8x^2-3` coefficient for coefficient |
| `14_43` | `INTSOL=1` | from the OBSTRUCTED class; produced under the integral-solution path. ⚠ The flag is recorded from a `ps` capture of the launch wrapper, not from the run log (lovelace's `genmodels.m` predates the line that prints it) — best available record, not log-confirmed | `ModelChecks` only (32 checks) — **no Guo-Yang equation exists for this base**, so there is no external oracle |

So, e.g.:

    CMNONCOPRIME=1 NORMALIZ_BIN=... magma -b D_s:=39 N_s:=2 OUTDIR:=... genmodels.m < /dev/null

## ⚠ The two reasons are NOT the same, and must not be conflated

* **`22_3`, `15_2`, `22_5` — ACCIDENTAL drift.** A guard was added *after* these files were
  committed. `Y2TWIST=1` decides the quadratic twist by Eichler-Selberg point count instead of
  dropping the cover, and restores them. This is a defect to repair: the right long-term fix is to
  make that selection the default once it is trusted.
* **`39_2`, `14_3`, `26_3` — DELIBERATE.** These exist only under `CMNONCOPRIME=1`, which is **off
  by default because it has no theoretical guarantee**: the `p | gcd(d,N)` local factor has no live
  implementation (`kappaminuszero` is dead code), and at `26_3` two non-coprime discriminants give
  provably wrong values.
  ⚠ **BUT THE STATED JUSTIFICATION FOR THE GUARD IS NOW IN DOUBT (2026-09-06).** `26_3` is the base
  whose two bad discriminants (`-267`, `-708`) are the whole reason the filter exists — and with
  the filter OFF, so those rows are admitted, `26_3` produces a model whose every quotient is
  isomorphic to Guo-Yang's, with the conic matching coefficient for coefficient. The two wrong
  values are evidently not load-bearing for the covers. Three bases now produce GY-matching models
  with the guard off (`39_2`, `14_3`, `26_3`) and none is known to be harmed by it.
  ⇒ A sweep of the 25 `X0_D_N.m` re-derivation tests under `CMNONCOPRIME=1` is the evidence needed
  to decide whether the guard should become off-by-default. Until that finishes, the flag stays
  opt-in. Their justification is the **published equation**, not regeneration.
  **Do not turn the flag on globally to make them "reproducible"** — that trades a documented gap
  for an undocumented risk on every base.

⇒ The target is **not** "everything regenerates by default". It is "every non-reproducing file has
a recorded reason and an independent validation". That is what this table and the model headers
encode.

## What checks what

| test | checks | runs |
|---|---|---|
| `tests/ModelChecks.m` | STORED models structurally — genus, Weil divisibility, Eichler-Selberg point counts. Independent of the Borcherds/Schofer path that produced them | CI, 85 files, 8889 checks |
| `tests/GuoYangEquations.m` | STORED models against the published equations, 10 bases | CI, ~125 s |
| `tests/_offline/ModelRegen.m` | that models still REGENERATE — the only check that runs the pipeline over stored files | offline |
| `tests/CRV_15_4.m` | `15_4`'s FULL genus-5 curve against trace-formula point counts — the only check of a `CRV` entry anywhere, and the one that pins its conic constant | CI, ~1 s |
| `tests/_offline/GuoYangCurve_14_3.m` | `14_3`'s full curve against Guo-Yang | offline, ~2 h |
| `tests/X0_D_N.m` (27 files) | re-derive the curve via `AllEquationsAboveCovers` and compare to stored/hand-written data — passing IS reproduction | CI |
| `tests/_offline/X0_87_1.m`, `X0_57_1.m`, `X0_14_5.m` | the same, for bases too slow for CI (`87_1` runs 45+ min) — `run_tests.m` globs only `tests/*.m`, so `_offline` is excluded automatically | offline |

`ModelRegen`'s `MR_KNOWN_DRIFT` lists exactly the six flagged bases above (the five originals
plus `14_43`).

⚠ **A KNOWN WEAKNESS OF THE `X0_D_N.m` TESTS.** `test_AllEquationsAboveCoversSingleCurve` SILENTLY
SKIPS cover keys it does not find (`if not is_def then continue`), so a base whose re-derived `W`
keys do not match the expected ones would pass **vacuously**. `X0_51_1` was negative-controlled by
hand (perturbing one coefficient makes it fail); the others have not been. The durable fix is a
comparison counter inside the helper that errors on zero — it would protect all 27 at once. **Update both it and this table
together**, or the next session gets a false "drifted" and repeats a day of work.

## ⚠ TWO MODELS WITH NO EXTERNAL ORACLE (`10_61`, `14_43`) — a weaker class of evidence

Added 2026-09-06 from the OBSTRUCTED class. Neither `X_0^10(61)` nor `X_0^14(43)` has a published
Guo-Yang equation, so **nothing outside our own pipeline confirms them**. `ModelChecks` alone
stands behind them (63 and 32 checks). Its strongest component is the Eichler-Selberg point count,
which *is* independent of the Borcherds/Schofer path that produced the models — but that is not the
same as matching a published curve. **Quote these two at that lower confidence**, and do not cite
them as evidence the pipeline is correct on obstructed bases in the way the Guo-Yang bases are.

## ⚠ PROCESS HAZARD: never `git pull` a clone that has jobs running from it

`10_61` and `14_43` were both started from lovelace's clone on 2026-09-04; that clone was pulled up
to `main` on 2026-09-06 **while both jobs were still running from it**. `AttachSpec` loads packages
**on demand**, so a long run can compile source that changed underneath it, and the output cannot
be pinned to one commit. No harm is evident (all checks pass), but reproducibility is the point of
this file. **Launch long runs from a COPIED tree** (as `~/shimura/vxfix` does) and leave the clone
free to update.

## ⚠ A FOURTH CATEGORY: literature-derived, NOT pipeline-produced (`15_4`)

`models_15_4.m` was **not** produced by this pipeline and never can be. Guo-Yang's published
Remark 39 (Compositio 153 (2017); absent from arXiv v1) says `X_0^15(4)`'s equation "is not
obtained using our method... the normalizer of the Eichler order in this case is larger than the
Atkin-Lehner group", and they take it from Tu (Pacific J. Math. 269 (2014), Lemma 13). Since
`N^+_B(O)` strictly contains `W_{15,4}`, the star quotient we form is the wrong object.

What we did instead: transcribe Tu's CM values (as quoted by GY) and **derive and check the rest**.
The polynomials are forced (they are the minimal polynomials of those values), our genus formula
independently reproduces GY's quotient structure, and the y-side constant `a = -1` is **confirmed
by our Eichler-Selberg point counts, discriminatingly** — 16/16 at `a = -1`, while `a = 1, 2, -2,
5, -5` each fail 2-3.

The conic constant `b = -1` is **also confirmed by us** (closed 2026-09-06; it was open for one
commit). It could not be reached the obvious way — `VerifyModelSet` is **vacuous** for `b`, passing
16/16 for *every* value, since the conic is genus 0 and point counts do not constrain its twist.
What closed it was the **full genus-5 curve**, which depends on `b` and which `ModelChecks` never
tests because it **skips `CRV` entries**: counting points on the fibre product over `F_p` and
comparing with `ComputePointsViaTrace` gives a match at all 12 primes 7..47 for `b = -1`, while
`-2,-3,-5,-6,-7,-10,-15,-30,1,2,3,5,15` each fail at 4-6. See `tests/CRV_15_4.m`.

⚠ **That exposed a general gap: no `CRV` (paired) entry has ever been point-count checked** — not
`93_1`, `57_1`, `21_2` either. Those bases are validated only through their hyperelliptic
quotients. Generalising `CRV_15_4.m` needs the ambient weights recorded in the model files, the
same gap noted above for `ModelRegen`.

⇒ `15_4` does **not** count toward "bases our pipeline reproduces". The Guo-Yang denominator for
that statistic stays **42**.

## A model that needs no flag but DOES need a code fix: `93_1`

`models_93_1.m` regenerates with the default recipe — but only on `main` **at or after the vx fix**
(`n_oo`, `BorcherdsForms.m` ~`:771`/`:787`/`:858`). Before it, the base died in the odd-`D`
`oo`-expansion block, which is why `93_1` had no model until 2026-09-05. Cost: 50927 s (14.1 h).

This is a **third** category, distinct from the two below: not a flag, not drift, but a *minimum
code version*. Anyone bisecting `main` to an older commit will find this file unregenerable and
should not read that as corruption. The same fix is what `95_1`, `115_1`, `123_1`, `129_1`,
`159_1` were blocked on.

## Reproducibility status, measured 2026-09-05

Of the 36 Guo-Yang bases we reproduce (`93_1` added 2026-09-05; the denominator is 42, not
43 — `15_4` is outside the method by their Remark 39):
* **24** are verified BY re-derivation (the `X0_D_N.m` tests run the pipeline, so passing IS
  reproduction);
* **11** are stored-model comparisons; of those, `51_1 57_1 14_5 55_1 21_2 87_1` regenerate clean
  (`ModelRegen`, measured), `93_1` regenerates clean given the vx fix (see above), and
  `39_2 14_3 22_3 15_2` need the flags above.
⚠ `93_1` is checked at the level of its `V_4` **quotients**, not its genus-5 full curve — the only
base in `GuoYangEquations.m` for which that is true. Its three cover keys `[1,3] [1,31] [1,93]` are
pinned against Guo-Yang; the `W={1}` entry is their fibre product and is validated only
structurally (`ModelChecks`). Closing that gap needs an offline full-curve test.
⚠ `21_2` and `57_1` report `OK` with **1 CRV skipped** — `ModelRegen` cannot rebuild `CRV` entries,
and for those two the `CRV` entry IS the `W={1}` full curve. So their *hyperelliptic covers* are
known to regenerate; their full curves are validated against Guo-Yang but not shown to regenerate.

✅ **The stated blocker for that is GONE (2026-09-06): the weights do NOT need recording.** They are
derivable — `y`'s weight is half the degree of its own equation, every other variable has weight 1 —
and `tests/CRVStructure.m` verifies that derivation reconstructs **16 of 21** stored `CRV` entries
exactly (irreducible, and the genus matches the one recorded beside it). No data migration is
needed; `ModelRegen` and the `X0_*` helper can derive the ambient space when they need it.

⚠⚠ **AND IT FOUND A REAL DEFECT: 5 `CRV` entries are DEGENERATE as stored.** `10_3` (four keys) and
`22_3` `[1,3]` each store the SAME equation twice, up to `y` <-> `x` — e.g.
`y^2 + 7/20*s^2 - 43/20*s*z + 2*z^2` alongside `x^2 + (the identical form)`. Then `y^2 = x^2`, the
scheme is REDUCIBLE (measured), and it cannot be the genus-1 curve recorded beside it. A fibre
product of a double cover with itself is reducible by construction, so this looks like two covers
with the same equation being paired as though independent.
**Nothing had ever checked these**, because `VerifyModelSet` skips every `CRV` entry.

✅ **ROOT-CAUSED AND FIXED the same day.** `EquationsAbovePointlessConics` builds the fibre product
`[y^2 - eqn2, x^2 - eqn1]`, taking `eqn2` from a cover of hyperelliptic degree `g+1` and `eqn1` from
a conic. **At `g = 1` the required degree `g+1 = 2` is also a conic's degree**, so the conic passed
the degree test and could be selected for BOTH roles — each of the five stored its own parent conic
twice (`10_3 [1,10]` against `[1,2,5,10]`, etc.). Fixed by requiring the two roles be filled by
different covers; those covers now DEFER, which is the honest outcome, and the five stored entries
were emptied to match. `CRV_KNOWN_BAD` is now EMPTY and must stay empty — a new degenerate entry is
a regression, not something to add to the list.
⚠ Cost of the repair: `10_3` loses 4 entries and `22_3` loses its `[1,3]` quotient. Neither costs an
external oracle — `10_3` is not a Guo-Yang base, and `22_3`'s Guo-Yang-validated `W={1}` entry is a
genus-3 hyperelliptic, untouched.
