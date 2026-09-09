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

## ⚠⚠ THE COPRIME-TO-LEVEL CM FILTER IS NOW **OFF BY DEFAULT** (flipped 2026-09-07)

`CMNONCOPRIME` is retired — what it used to enable is the default. **`CMCOPRIME=1` restores the old
filtering** if a base ever turns out to be poisoned by an admitted point.

**Evidence for the flip.** A full sweep of the 11 `N>1` `X0_D_N.m` re-derivation tests, each run
both ways: **10 of 10 pass identically**. (For `N = 1` the filter is provably a no-op, excluding 19
of 30 tests rigorously rather than by sampling.) Two tests appeared to fail with it off — `10_13`
and `6_17` — and both were artifacts of a hardcoded coordinate matrix *in the test*, gone once the
isomorphism is constructed rather than pinned. Decisively, `26_3` is the very base whose two
misbehaving discriminants justified the filter, and with them admitted its full `V_4` diagram still
matches Guo-Yang, conic coefficient for coefficient. The filter also **cost** models: at `bd := 4`
it cut `26_3`'s pool from 21 to 3 against demand 15, and `39_2`'s from 24 to 3 against 19.

**Measured consequences** (verified by regeneration, not assumed):

| base | before | after |
|---|---|---|
| `26_3` | needed the flag | builds with **no flag**, 15 keys |
| `14_3` | needed the flag | **no flag**, 15 keys, **identical to committed** |
| `51_1` | never needed it | still **byte-for-byte identical** |
| `39_2` | needed the flag | **no flag**, 15 keys, **identical to committed** |

⚠ **THE GAP THIS LEAVES OPEN.** There is **no theoretical guarantee**, only the sweep above. The
local factor at `p | gcd(d,N)` has **no live implementation** (`kappaminuszero` is dead code), and
Schofer's Thm 4.1 assumes unimodularity at unramified primes, which fails at a level prime where the
order is Eichler. The two known-wrong values at `26_3` (`-267`, `-708`) **are still wrong** — they
simply do not propagate into the cover equations. Therefore:
* a model built from non-coprime discriminants still needs an **INDEPENDENT ORACLE** (a published
  equation, or Eichler-Selberg point counts); regeneration alone is not enough;
* **do not read this flip as evidence the `p | gcd(d,N)` factor is unnecessary.** Supplying it is
  still the real fix — it is what would make the swap class *correct* rather than merely *harmless*.

## Files needing a non-default flag

| base | flag | why | validated by |
|---|---|---|---|
| `26_3` | *(no flag)* but `base_label := 8103` | stores the presentation whose `V_4` is the one Guo-Yang use, so the full curve is directly comparable. A DEFAULT run yields a different (equally valid) `V_4` and so differs in the `W={1}` entry only | `tests/CRVFullCurve.m` — full-curve isomorphism CONSTRUCTED and certified |
| `14_43` | `INTSOL=1` | from the OBSTRUCTED class; produced under the integral-solution path. ⚠ The flag is recorded from a `ps` capture of the launch wrapper, not from the run log (lovelace's `genmodels.m` predates the line that prints it) — best available record, not log-confirmed | `ModelChecks` only (32 checks) — **no Guo-Yang equation exists for this base**, so there is no external oracle |

Everything else uses the plain recipe above. (`CMNONCOPRIME` no longer exists as a flag — what it
enabled is the default; `CMCOPRIME=1` is the escape hatch in the other direction.)

## ⚠ The two reasons are NOT the same, and must not be conflated

* **`22_3`, `15_2`, `22_5` — RESOLVED 2026-09-07, and NOT the way this file predicted.** They no
  longer need any flag: regenerated with the plain recipe they give MORE covers than the files they
  replaced (`22_5` 3 → 11 populated, `15_2` 12 → 15, `22_3` 13 → 15), nothing lost, and
  `GuoYangEquations` still passes.
  ⚠ **The fix was the COPRIME flip, not `Y2TWIST`.** This file used to say "the right long-term fix
  is to make that selection the default once it is trusted". `Y2TWIST` was evaluated for exactly
  that and left OFF: a controlled run — default vs the selector disabled, on the SAME code — is
  IDENTICAL at all three bases, and the deferral path logs zero "unpinned y2-scale" messages. It
  never fires here any more.
  ⚠⚠ The first evaluation got this backwards by comparing `Y2TWIST=1` runs against the COMMITTED
  files, which predate the coprime flip — so the coprime flip's gains were credited to `Y2TWIST`.
  **Compare against a current baseline, never a committed artifact.**
* **`39_2`, `14_3`, `26_3` — NO LONGER NEED A FLAG (2026-09-07).** The coprime filter is now off by
  default, and all three regenerate without one (`39_2` and `14_3` byte-identical to committed).
  What follows is kept because the underlying THEORETICAL gap is unchanged: the `p | gcd(d,N)` local factor has no live
  implementation (`kappaminuszero` is dead code), and at `26_3` two non-coprime discriminants give
  provably wrong values.
  ⚠ **THAT DOUBT IS RESOLVED — the sweep ran and the guard was flipped (see the section above).**
  `26_3` is the base whose two bad discriminants (`-267`, `-708`) were the filter's whole
  justification, and with them admitted its every quotient is isomorphic to Guo-Yang's, the conic
  coefficient for coefficient. Their justification remains the **published equation**, not
  regeneration.

⇒ The target is **not** "everything regenerates by default". It is "every non-reproducing file has
a recorded reason and an independent validation". That is what this table and the model headers
encode.

## What checks what

| test | checks | runs |
|---|---|---|
| `tests/ModelChecks.m` | STORED models structurally — genus, Weil divisibility, Eichler-Selberg point counts. Independent of the Borcherds/Schofer path that produced them | CI, 88 files, 9349 checks |
| `tests/GuoYangEquations.m` | STORED models against the published equations, 11 bases | CI, ~122 s |
| `tests/_offline/ModelRegen.m` | that models still REGENERATE — the only check that runs the pipeline over stored files | offline |
| `tests/CRVFullCurve.m` | CRV pairs against Guo-Yang by CONSTRUCTED full-curve isomorphism — Mobius map from the hyperelliptic quotient, then `IsIsomorphism` certifies it. Proof, not a screen; avoids the generic call that runs for hours on these | CI, ~0.1 s |
| `tests/CRV_15_4.m` | `15_4`'s FULL genus-5 curve against trace-formula point counts — the only check of a `CRV` entry anywhere, and the one that pins its conic constant | CI, ~1 s |
| `tests/_offline/GuoYangCurve_14_3.m` | `14_3`'s full curve against Guo-Yang | offline, ~2 h |
| `tests/GuoYangQuotientOracle.m` | **144** quotients over **20 bases**, derived generically via `CurveQuotient` from Guo-Yang's published curve + involutions | CI, ~4.6 s |
| `tests/GuoYangQuotients_10_19.m` | **12** `10_19` quotients, hand-derived (CRV pair — `CurveQuotient` is blocked, see below) | CI, 0.26 s |
| `tests/GuoYangQuotients_22_5.m` | **9** `22_5` quotients, hand-derived (non-diagonal action) | CI, 0.13 s |
| `tests/GuoYangQuotients_10_13.m` | **9** `10_13` quotients, hand-derived; determines a labelling error in Guo-Yang's table | CI, 0.05 s |
| `tests/_offline/FullCurve_22_5.m` | CONSTRUCTS `22_5`'s full curve; reproduces Guo-Yang verbatim | offline, ~420 s |
| `tests/X0_D_N.m` (34 files) | re-derive the curve via `AllEquationsAboveCovers` and compare to stored/hand-written data — passing IS reproduction | CI |
| `tests/_offline/X0_87_1.m`, `X0_57_1.m`, `X0_14_5.m` | the same, for bases too slow for CI (`87_1` runs 45+ min) — `run_tests.m` globs only `tests/*.m`, so `_offline` is excluded automatically | offline |

## THE QUOTIENT ORACLE — the largest source of external validation here

⚠ **`tests/GuoYangEquations.m` compares only the equations Guo-Yang PRINT**, which for most bases
is the full curve alone — so a base with fifteen cover keys got **one** external comparison. But
they also print the **involutions**, and every quotient follows from those:
`CurveQuotient(AutomorphismGroup(C,[w]))` is `X/W`. That turns one comparison per base into one per
cover key, against a source entirely outside the pipeline. Current totals: **174 quotient
comparisons over 23 bases**, all passing, in about 5 seconds.

**Three errors in Guo-Yang's tables have now been determined this way or alongside it**, each by
evidence rather than by preferring one source:
* `93_1`: `-3t` is a typo for `-3s` — confirmed independently by the journal version.
* `14_5`: the table's `w_35` sign is wrong; **their own Example 36** has it right, and the table's
  map is `w_10`.
* `10_13`: the table **swaps `w_10` and `w_13`** — settled by Ogg's fixed-point rule (the fixed
  points of `w_m` are the CM points of discriminant `-4m`), and the clincher is internal: **their
  own CM table** puts disc `-52` at Hauptmodul `0` and `-40` at infinity, contradicting their
  involution table and agreeing with our pipeline.

⚠ **WHAT THE ORACLE CANNOT COVER, and why.** `CurveQuotient` fails on any curve whose ambient is a
weighted projective space — Magma models that as a toric variety and `IdentityMap` returns a
`TorMap`, so `AutomorphismGroup` dies (reported as
[Magma-Maths/Magma#123](https://github.com/Magma-Maths/Magma/issues/123)). That is exactly the
`CRV` paired presentations, so `10_19`, `22_5` and `10_13` each need a hand-derived oracle file.
`10_23` is excluded for a different reason: genus 9, where `CurveQuotient` ran >11 min and got the
process OOM-killed.

⚠ **TRANSCRIPTION TRAPS, all of which bit.** `pdftotext` drops superscripts, so `10_23`'s `w_2`
reads `-55 y` for `-5^5 y`, and every `y/x6` is `y/x^6`. And **the automorphism check cannot catch
a dropped sign when `f` is even**: at `6_29` I wrote `w_3 = (2/x, 8y/x^6)` for their
`(-2/x, 8y/x^6)`, and since `f(2/x) = f(-2/x)` both preserve the curve — but they differ by `w_2`,
so the quotients came out swapped. The **per-label genus comparison** is what caught it.

`ModelRegen`'s `MR_KNOWN_DRIFT` lists exactly the six flagged bases above (the five originals
plus `14_43`).

✅ **CLOSED: ALL 34 `X0_*.m` TESTS NOW CHECK INVOLUTIONS** (was 23 on 2026-09-07).
The tests generated on 2026-09-07 all carried an EMPTY `ws_data`, so they made zero involution
comparisons: they verified each cover is isomorphic to the stored curve, but not that the
Atkin-Lehner involutions correspond — and the involutions are what make these QUOTIENT models
rather than merely curves.

**All of `51_1 55_1 22_3 15_2 14_5 26_3 57_1 21_2 14_3 15_1 22_5`** now check involutions, so the
count is **34 of 34**. `22_5` was last and needed the pipeline change below: before
`EquationsByRebase` its `W={1}` cover did not exist at all, so there was nothing to attach
involutions to. Its expected curve is Guo-Yang's own published equation, so their involutions apply
verbatim with no transport.
For `26_3` and `57_1` this also added the `W={1}` CRV pair itself to `cover_data`, which the
generator had omitted; `psi` there comes from `construct_crv_isomorphism` rather than
`IsIsomorphic`, which hangs on paired presentations. The matrices
were obtained by TRANSPORT, which is what makes them non-circular: Guo-Yang publish the
involutions in THEIR coordinates, `psi := IsIsomorphic(our stored curve, their curve)` is computed
from the two EQUATIONS alone — never from the pipeline's own `ws` — and the recorded matrix is
`psi^-1 . w_GY . psi`, which came out linear in the weighted coordinates in all 11 cases. The
script is `tests/_gyinvol.m`; each matrix is checked to be an involution of our curve and to equal
the transported map. `psi` is one element of a torsor under `Aut`, but the harness searches that
same torsor, so the choice cannot produce a false verdict either way. **Negative-controlled**: on
`51_1`, swapping `w_3` and `w_51` makes the test fail on the labelling — and on `14_5` this
machinery *determined a typo in the journal's table* (see below).

⚠ **`26_3`'s test now REQUIRES `base_label := 8103`**, and this is not cosmetic. `models_26_3.m`
deliberately stores the `V_4` Guo-Yang use, which a default run does not produce — it gives a
different, equally valid one — so without the label the `W={1}` pair the pipeline emits is a
genuinely different presentation and the isomorphism assertion fails.

⚠ **`14_3` NEEDED `base_label := 5394`, AND FINDING IT IS THE LESSON.** The default base (5383)
gives a valid pair that is NOT Guo-Yang's `V_4`; `construct_crv_isomorphism` declines, and the
general `IsIsomorphic` needs **6739 s (112 min)** to confirm the curves are abstractly isomorphic
while still producing no usable coordinate change. Sweeping the candidate bases — replaying only
the pointless-conic step, ~7 s each — finds 5394, whose pair IS their `V_4`:

|  | `y`-equation | conic |
|---|---|---|
| ours (base 5394) | `y^2 = 4s^4 + 88s^2z^2 - 28z^4` | `x^2 = -2s^2 - 9z^2` |
| Guo-Yang | `y^2 = -7x^4 + 22x^2 + 1` | `z^2 = -9x^2 - 2` |

theirs scaled by 4 with the two base coordinates exchanged. The constructor then succeeds
immediately and the linear solve drops from 70 monomial equations to 10.
⇒ **SWEEP THE BASE BEFORE CONCLUDING ANYTHING ABOUT A PAIR** (same lesson as `26_3`'s 8103).
Tool: `tests/_basesweep.m`.

⚠ **`21_2` IS DIFFERENT, AND THE SWEEP DOES NOT RESCUE IT.** All four candidate bases were tried
and none reproduces Guo-Yang's `V_4`, because theirs has `y` of weight 3 (a genus-2 `y`-quotient)
while every base our pipeline offers gives weight 2 (genus 1). So for `21_2` the slow
`IsIsomorphic` + linear-solve route was genuinely necessary, not merely a route not yet optimised.

✅ **`22_5` IS NOW CLOSED**, and the reason it was open is worth keeping.

⚠ **WHAT GUO-YANG DO DIFFERENTLY AT `22_5`: THEY DO NOT BUILD IT AS A FIBRE PRODUCT.** Their model
is a single plain hyperelliptic curve, `y^2 = -11x^12 - 80x^10 - 240x^8 - 362x^6 - 240x^4 - 80x^2
- 11` (degree 12, genus 5), a double cover of `X/w_110` with `w_110 = (x,-y)` hyperelliptic. The
paper gives `22_5` no worked example and flags no exception for it (the only one anywhere is
Remark 39, `15_4`), so it came out of their standard method.

Our pipeline never constructs a genus-`g` curve directly. `process_P1_cover` assembles it as a
FIBRE PRODUCT, needing another quotient to carry an equation of degree exactly `g+1` **over a base
shared with the `P^1`/conic**. At `22_5` that fails twice and cascades: the genus-2 quotients
`{1,2}`, `{1,5}`, `{1,11}` each need **degree 3** and the degrees produced are `1,2,4,6,7,8`; and
`{1}` needs **degree 6** over a base shared with the `P^1` `[1,110]`, while the only degree-6
equation sits over the star base.

✅ **SOLVED 2026-09-08, AND IT IS A ONE-LINE IDEA: CHANGE THE HAUPTMODUL ON THE STAR BASE.** Which
degree a quotient's equation has depends on whether the point at infinity is a branch point — i.e.
on the Hauptmodul normalisation, which is ours to choose. The Mobius change `t -> r + 1/u`, for `r`
a **rational root** of an equation over that base, moves a branch point to infinity: a quartic
becomes a cubic, and a conic that shares the root collapses to degree 1, a `P^1`. At **`r = 4/5`**
both happen at once, for `{1,2,11,22}` (4 -> 3) and `{1,10,11,110}` (2 -> 1) — precisely the two
covers `{1,11}` needs. Once the genus-2 quotients exist, `{1}` follows.

**The result is Guo-Yang's polynomial VERBATIM, coefficient for coefficient**, not merely
isomorphic to it — and nothing in the construction knows their equation. `r = 0` builds a model
that is isomorphic but NOT equal, so the coincidence at `r = 4/5` is real rather than an artefact
of any construction landing on genus 5. `tests/FullCurve_22_5.m` pins both (~420 s), and also
recovers all three genus-2 quotients.

⚠ **BUT IT DOES NOT GENERALISE, AND THAT WAS TESTED RATHER THAN ASSUMED.** At `10_19` the same
sweep fills 3 of 4 empty keys, yet the genus-5 pair it produces has a `y`-side that is **not
isomorphic to any Atkin-Lehner quotient of `X_0^10(19)`** — checked against all three genus-2
quotients derived from Guo-Yang's own curve and involutions. Its conic IS right, which is precisely
why the output looks plausible; **filling a key is not building the right curve**. So the lever is
correct at `22_5` (verified verbatim against a published equation) and WRONG at `10_19`, and no
output of `tests/_rebaselever.m` should be accepted without an independent oracle. Why it diverges
is not diagnosed; the likely suspect is that the rebase changes the coordinate while the `y^2`-scale
(twist) is not re-derived from CM values afterwards.

NOT wired into `AllEquationsAboveCovers` — and on this evidence it should not be until the failure
mode is understood. `models_22_5.m` still records the empty entries and `tests/X0_22_5.m` still
checks no involutions.

✅ **A REAL WIN CAME OUT OF TESTING IT, THOUGH: A COMPLETE ORACLE FOR `10_19`.** Guo-Yang print only
three quotients for that base, so only three of our fifteen keys looked checkable — but they also
print the full curve AND the involutions, and every quotient follows from those. The action is
diagonal, so each quotient is the field of invariant monomials. `tests/GuoYangQuotients_10_19.m`
checks **11 of our stored entries** against it, label by label, in 0.24 s: **11 matches, 0
mismatches**. The construction is calibrated against the two quotients Guo-Yang state in words
(`X/w_190` genus 2, `X/w_38` genus 0) before being trusted. ⇒ **This technique applies to every base
where they print the top curve with its involutions**, and is much cheaper than re-deriving
anything. `22_5` has an EMPTY `[1]` entry in its model file, so there is
no full curve to attach Guo-Yang's involutions to at all — that one is structural, not effort.
**Counted, not assumed** (`SetVerbose("ShimuraQuotients",1)` prints them; the repo has produced
three vacuous tests, so the comparisons made are checked rather than inferred from a green run):

| base | curve cmps | involution cmps | covers matched |
|---|---|---|---|
| `51_1` | 4 | 2 | 4/4 |
| `55_1` | 4 | 2 | 4/4 |
| `22_3` | 14 | 3 | 14/14 |
| `15_2` | 13 | 3 | 13/13 |
| `14_5` | 8 | 4 | 8/8 |
| `26_3` | 12 | 3 | 12/12 |
| `57_1` | 4 | 2 | 4/4 |

19 involution comparisons in total, and every expected cover matched in every case.

`14_3` and `21_2` needed a different route, and it is now the general one. `construct_crv_isomorphism`
DECLINES on both, because our pair and Guo-Yang's present the curve over DIFFERENT intermediate
quotients (at `21_2` Guo-Yang's `y` has weight 3 and a genus-2 `y`-quotient, ours weight 2 and
genus 1), so there is no common base to take a Mobius map from. The general `IsIsomorphic` is the
fallback: 208 s at `21_2`, **6739 s (112 min) at `14_3`**, both returning true.

⚠ **THREE TRAPS ON THIS ROUTE, each of which cost a wrong conclusion or a rerun:**
1. `Inverse(psi)` raises `"Map has no inverse"` on the very map `IsIsomorphic` returns.
   `IsInvertible` succeeds on the SAME map — a representation issue, not a mathematical one.
2. The composite `psi^-1 . w_GY . psi` comes back as ONE unreduced degree-39 representation
   (`AllDefiningPolynomials` offers no other), so **reading coefficients off it reports "not
   linear" although the MAP is linear**. Recording that as non-linearity would have been a
   wrong-object claim of exactly the kind this repo keeps paying for.
3. The fix is to **SOLVE for the matrix, not read it**. On `P(1,wy,1,1)` only `y` has weight `wy`,
   so a weight-respecting matrix must send `y -> c*y` and act on `(x,s,z)` by a 3x3 block; and
   `q1*L3-q3*L1` and `q1*L4-q4*L1` vanishing on the curve are LINEAR conditions on that block's 9
   coefficients. At `21_2` the kernel came out **dimension 1** for all three involutions — the
   block is unique up to scalar — and each was certified by MAP EQUALITY against the transported
   map, which is representation-independent. Script: `tests/_gyinvol_crv.m`.

⚠ **A SECOND GUO-YANG TABLE TYPO DETERMINED, at `14_5`.** The journal's table prints
`w_35(x,y) = ((x+2)/(2x-1), -25y/(2x-1)^4)` while its own Example 36 prints `+25y`. Both are
involutions of the curve, so inspection cannot choose between them; they differ by `w_14`, and
`35*14/gcd(35,14)^2 = 10`, so the two readings are `w_35` and `w_10`. `tests/X0_14_5.m` adjudicates:
`+25` labelled `w_35` passes, `-25` labelled `w_35` FAILS on the labelling, `-25` labelled `w_10`
passes. Example 36 is right and the table is wrong. Both readings are kept in the test under their
own labels, so the run makes 4 involution comparisons and neither can be quietly relabelled.

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

Of the 38 Guo-Yang bases we reproduce (`93_1` added 2026-09-05; the denominator is 42, not
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
