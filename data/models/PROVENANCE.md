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

So, e.g.:

    CMNONCOPRIME=1 NORMALIZ_BIN=... magma -b D_s:=39 N_s:=2 OUTDIR:=... genmodels.m < /dev/null

## ⚠ The two reasons are NOT the same, and must not be conflated

* **`22_3`, `15_2`, `22_5` — ACCIDENTAL drift.** A guard was added *after* these files were
  committed. `Y2TWIST=1` decides the quadratic twist by Eichler-Selberg point count instead of
  dropping the cover, and restores them. This is a defect to repair: the right long-term fix is to
  make that selection the default once it is trusted.
* **`39_2`, `14_3` — DELIBERATE.** These exist only under `CMNONCOPRIME=1`, which is **off by
  default because it has no theoretical guarantee**: the `p | gcd(d,N)` local factor has no live
  implementation (`kappaminuszero` is dead code), and at `26_3` two non-coprime discriminants give
  provably wrong values. Their justification is the **published equation**, not regeneration.
  **Do not turn the flag on globally to make them "reproducible"** — that trades a documented gap
  for an undocumented risk on every base.

⇒ The target is **not** "everything regenerates by default". It is "every non-reproducing file has
a recorded reason and an independent validation". That is what this table and the model headers
encode.

## What checks what

| test | checks | runs |
|---|---|---|
| `tests/ModelChecks.m` | STORED models structurally — genus, Weil divisibility, Eichler-Selberg point counts. Independent of the Borcherds/Schofer path that produced them | CI, 82 files, 8767 checks |
| `tests/GuoYangEquations.m` | STORED models against the published equations, 9 bases | CI, ~97 s |
| `tests/_offline/ModelRegen.m` | that models still REGENERATE — the only check that runs the pipeline over stored files | offline |
| `tests/_offline/GuoYangCurve_14_3.m` | `14_3`'s full curve against Guo-Yang | offline, ~2 h |
| `tests/X0_D_N.m` (25 files) | re-derive the curve via `AllEquationsAboveCovers` and compare to hand-written data | CI |

`ModelRegen`'s `MR_KNOWN_DRIFT` lists exactly the five bases above. **Update both it and this table
together**, or the next session gets a false "drifted" and repeats a day of work.

## Reproducibility status, measured 2026-09-05

Of the 34 Guo-Yang bases we reproduce:
* **24** are verified BY re-derivation (the `X0_D_N.m` tests run the pipeline, so passing IS
  reproduction);
* **10** are stored-model comparisons; of those, `51_1 57_1 14_5 55_1 21_2 87_1` regenerate clean
  (`ModelRegen`, measured), and `39_2 14_3 22_3 15_2` need the flags above.
⚠ `21_2` and `57_1` report `OK` with **1 CRV skipped** — `ModelRegen` cannot rebuild `CRV` entries
because the model file does not record the ambient weights, and for those two the `CRV` entry IS
the `W={1}` full curve. So their *hyperelliptic covers* are known to regenerate; their full curves
are validated against Guo-Yang but not shown to regenerate. Closing that gap means recording the
weights (recoverable: `y`'s weight is half the degree of its defining polynomial).
