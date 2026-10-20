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
suspect discriminants justified the filter, and with them admitted its full `V_4` diagram still
matches Guo-Yang, conic coefficient for coefficient. The filter also **cost** models: at `bd := 4`
it cut `26_3`'s pool from 21 to 3 against demand 15, and `39_2`'s from 24 to 3 against 19.

**Correction (2026-09-26): the `26_3` values at `-267`, `-708` appear to be right.** They were
recorded as an `s <-> s~` swap against Guo-Yang arXiv v1 Table 49; that v1 table appears instead to
have a sign misprint. Our `s(-267) = 8/25`, `s(-708) = 11/49` are `phi(9/25)`, `phi(27/49)` under
`phi(w) = (1-w)/2`, i.e. the sign-flipped printed values, while the printed `-9/25`, `-27/49` map to
`17/25`, `38/49`. Shimura reciprocity on Guo-Yang's own equations rejects the printed values and
accepts the flipped ones (`tests/CMPoints.m`, `(26,3)` block). The journal version (Compositio 153,
2017) has no CM-value tables, so no published table arbitrates.

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
order is Eichler. Therefore:
* a model built from non-coprime discriminants still needs an **INDEPENDENT ORACLE** (a published
  equation, or Eichler-Selberg point counts); regeneration alone is not enough;
* **do not read this flip as evidence the `p | gcd(d,N)` factor is unnecessary.** Supplying it is
  still the real fix — it is what would make the swap class *correct* rather than merely *harmless*.

## Files needing a non-default flag

| base | flag | why | validated by |
|---|---|---|---|
| `26_3` | *(no flag)* but `base_label := 8103` | stores the presentation whose `V_4` is the one Guo-Yang use, so the full curve is directly comparable. A DEFAULT run yields a different (equally valid) `V_4` and so differs in the `W={1}` entry only | `tests/CRVFullCurve.m` — full-curve isomorphism CONSTRUCTED and certified |
| `14_43` | `INTSOL=1` | from the OBSTRUCTED class; produced under the integral-solution path. ⚠ The flag is recorded from a `ps` capture of the launch wrapper, not from the run log (lovelace's `genmodels.m` predates the line that prints it) — best available record, not log-confirmed | `ModelChecks` only (32 checks) — **no Guo-Yang equation exists for this base**, so there is no external oracle |

| `21_1` | `HMFIT=1` | the hauptmodul normalisation is read off the two discriminants where a value VANISHES, so those satisfy the relation BY CONSTRUCTION and can never be flagged; at `21_1` the datum it reads FROM is the wrong one (five discriminants agree on `scale_tilde = 36`, the default takes `9` from `d = -7`). `HMFIT=1` solves for the pair against every rational CM point instead. ⚠ It also needed the fitted scale's SIGN pinned positive — see below | **`tests/GonzalezRotger.m`, an EXTERNAL oracle**: the `W=[1]` full curve has Jacobian `21a2`, matching the published equation (ours `-343x^4+94x^2-7`, theirs `-7x^4+94x^2-343` — the same curve by `x -> 1/x`). Negative-controlled: the `-1` twist gives `336e4` and the oracle rejects it. Also `VerifyModelSet` 44/0 (negative control: twisting the `[1,3]` entry gives 3 failures) and `ModelChecks` 44/0 |

**`33_1` and `69_1` (2026-09-14) use the PLAIN RECIPE — no flag.** Both were blocked by the
`ScaleForSchofer` bug at `d = -4` (see `HANDOFF.md`), not by anything about their data. `33_1` was
built twice, with and without `HMFIT=1`, and the two files are BYTE-IDENTICAL: it does not need the
flag. `33_1` is confirmed by `tests/GonzalezRotger.m` (Jacobian `33a1`, the published equation);
`69_1` has no external oracle and rests on `VerifyModelSet` 44/0 plus a negative control.

⚠ **`HMFIT` is still OFF by default and is validated at ONE base.** `21_1` is the first external
confirmation that its fit picks the right datum — real evidence, but a single base. It does NOT
unblock `33_1`, the other Gonzalez-Rotger target: that one still dies on the `Log11` runaway, a
different and still-open cause.

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
  implementation (`kappaminuszero` is dead code). (The `26_3` values at `-267`, `-708` once
  recorded here as provably wrong appear to be right; see the 2026-09-26 correction above.)
  ⚠ **THAT DOUBT IS RESOLVED — the sweep ran and the guard was flipped (see the section above).**
  `26_3` is the base whose two suspect discriminants (`-267`, `-708`) were the filter's whole
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
| `tests/_offline/X0_10_19.m`, `X0_21_1.m` | the same, for the two bases still too slow for CI (`87_1`, `111_1`, `93_1` and `39_2` moved into CI on 2026-10-01, minutes each with the Borcherds search of #63) — `run_tests.m` globs only `tests/*.m`, so `_offline` is excluded automatically | offline |

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

⚠ **MEASURED 2026-09-09: THE `X0_*` TESTS RE-DERIVE ONLY 41% OF THE COVERS.** Across the 34 tests
there are **128 `cover_data` keys against 309 populated model keys**. Many check a single cover:
`10_11 10_13 10_23 6_11 6_17 6_19 6_29 6_31 6_37` are 1 of 15 each, and nine level-one bases
(`134_1 194_1 35_1 38_1 39_1 58_1 62_1 74_1 82_1 86_1 94_1`) are 1 of 4. The helper SILENTLY SKIPS
a key that is absent from `cover_data` (`if not is_def then continue`), so this is invisible unless
the keys are counted.

⚠ **What this does and does not mean.** The stored MODELS are well checked — the quotient oracle
makes ~190 comparisons over 25 bases against Guo-Yang. What is thin is the RE-DERIVATION claim,
that the pipeline reproduces them: for most bases that rests on one cover. Closing it is mechanical
(emit `cover_data` from the model) but must handle two things that have already caused defects:
entries stored as `<genus, f, h>`, and `CRV` pairs. It also costs CI time — each key is an
`IsIsomorphic`.

**Done so far**: `206_1` went 1 -> 4 of 4 (its `[1,103]` entry carries `h = x^5+x^4+x^3+x^2`), and
`_offline/X0_87_1.m`'s long-standing failure was diagnosed as exactly this dropped `h`.

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

## `95_1`, `119_1`, `159_1`: default flags, but the Borcherds search of PR #63 (`zero-side-kernel`)

`models_95_1.m` was produced 2026-10-01 on lovelace with the default recipe (`BFPROGRESS=1
BFCACHE=1`, no result-changing flag) from the branch of PR #63 at `fa90a45`, the class-number
tables mounted (`CLASS_GROUPS_FAST_DIR=/scratch/class-groups-fast`). On `main` before that PR the
second rung of the Borcherds ladder (0-side pole order 3325) had never finished (7 h+); with it the
five forms take 23 min on a Mac and the ladder ends on the third rung (m = -43). The whole run,
all six stages, took about 2 h on a loaded lovelace. The top curve is isomorphic to Guo-Yang's
published equation (`tests/GuoYangEquations.m`), and genus 7 and the hyperelliptic involution match
Table 1. The quotient comparison (`tests/GuoYangQuotientOracle.m`) has no entry for `95_1`, so it checks
nothing here; the quotient models rest on the L-polynomial comparison against the trace formula in
`VerifyModelSet` and on the involution check in `tests/_offline/X0_95_1.m`.

`models_119_1.m` was produced the same way, the same day: its ladder ends on the third rung
(m = -51, 0-side pole order 6069, 52 good forms), the whole run took about 2 h 10 min on a loaded
lovelace, and the top curve (genus 9) is isomorphic to Guo-Yang's published equation and matches
Table 1. The quotient comparison has no entry for `119_1` either; its quotients are checked the same
way as `95_1`'s.

`models_159_1.m`, the same way and day: its ladder ends on the third rung (m = -51, 0-side pole
order 8109, pool 8234, 52 good forms), whose kernel alone took about 4 h 40 min on a loaded
lovelace (5 h 40 min for the whole run); the top curve (genus 9) is isomorphic to Guo-Yang's
published equation (`tests/GuoYangEquations.m`), genus and hyperelliptic involution match Table 1,
and `tests/_offline/X0_159_1.m` checks Guo-Yang's w_3. The quotient comparison has no entry for
`159_1`, so it checks nothing there. Its `[1,53]` quotient is stored with an `h` term
(`y^2 + h y = f`).

For all three, every stored quotient's L-polynomial at the small good primes equals the one the
trace formula gives for the W-fixed part of the D-new space (`tests/ModelChecks.m`, check [5]),
which also fixes the labels of the w_p: the opposite sign convention would swap the genera of the
two middle quotients at each base. The re-derivation tests for `95_1` and `119_1` take 28 and
58 min on a Mac and live in `tests/_offline/` with `159_1`'s.

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

## 203 covers built as fibre products of their double covers (2026-10-03/04, `FibreProductCovers.m`)

A cover `X_0(D,N)/W` is Galois over the star curve with an elementary abelian 2-group, so it is
the fibre product over the star line of its index-2 Atkin–Lehner double covers `y^2 + hy = f(t)`,
i.e. `(y + h/2)^2 = f + h^2/4`. The entries added on these dates, 203 of them across 53 files, were
built that way from the hyperelliptic entries already committed in the same file, with no Borcherds
or Schofer computation, by the writer `tools/writefp.m` (the same construction the pipeline's final
stage runs). Among them are 29 top curves `W = {1}`, of genus 4 to 11; the rest are intermediate
quotients of genus 2 to 11.

⚠ **At level `N > 1` these entries rest on the `m = 0` term of Schofer's formula**, through the
double covers they are built from, and PR #66 only *proposes* a proof of that term; so they are
checked by point counts against the trace formula, not proved. The same goes for the three
"not hyperelliptic" verdicts the construction gives (`6_23/w_23`, `34_3/w_2`, `46_3/w_3`).

**Format.** Stored as the existing non-hyperelliptic shape `<genus, "CRV", [ Strings() | ... ]>`,
one equation per factor, in coordinates `s, z` of weight 1 and fibre coordinates `y1, y2, ...`
each of weight half the degree of its own equation (the same rule the older `x, y, s, z` pairs
follow). Twenty-nine of them are top curves `W = {1}`, of genus 3 to 11, which no earlier stage
could reach because they are not subhyperelliptic.

**What vouches for them.** Each entry was accepted only when the compositum has the genus the
Shimura-curve genus formula predicts and its point counts over `F_p` and `F_{p^2}` at the first
three good primes agree with the Eichler–Selberg trace formula on the `W`-fixed part of the `D`-new
space. A prime at which a factor reduces badly is skipped, so eleven entries rest on two primes
rather than three: `14_19` `[1]`, `[1,7]`, `[1,14]`; `35_2` `[1]`, `[1,10]`, `[1,14]`; `62_3` `[1]`,
`[1,2]`, `[1,6]`, `[1,31]`, `[1,93]` (a prime where some factor has a repeated root mod `p`; with
the writer's weaker criterion, a factor that is a square mod `p`, only the three `62_3` entries
`[1]`, `[1,2]`, `[1,6]` do). Every written entry was then re-read from its stored strings alone and
re-checked the same way at two primes. Nothing external publishes these curves, so there is no
published-equation comparison; the two genus-2 curves at `34_3` are isomorphic to Saia's models
(Sachi, review of PR #72), and all 200 of the first batch agree with point counts from modular
symbols at 9–13 primes each (same review).

**Rejections, and the two counts.** The writer of 2026-10-03 tried one generating set of double
covers per cover and rejected 16 candidates (14 by genus, 2 by the trace formula). The stage's own
sweep (PR #71) tries every generating set and reported 11 failures over a pool that also held the
24 stored pairs: the difference is six covers that a second generating set builds. The writer now
tries every set too (and skips a degenerate one, where some subset of factors multiplies to a
constant times a square), which reaches seven more covers; three are added here, on 2026-10-04,
`22_7` `[1,7]`, `26_5` `[1,13]` and `74_3` `[1,111]`, and the other four have models from reruns
instead: `21_2` `[1,3]`, `[1,6]` (PR #70) and `26_5` `[1,10]` (PR #73) as `y^2 = f(x)`, and `22_7`
`[1,22]` (PR #73) as a genus-3 double cover of the conic `{1,14,22,77}`, which has no rational
point, so that curve has no model `y^2 = f(x)` over Q. What remains rejected, with
every generating set: `14_19` `[1,38]` (`[1,19]` has a hyperelliptic model from the rerun below); `34_7` `[1,7]`, `[1,119]`; `82_3` `[1,2]`,
`[1,82]`, `[1,123]` (no set has the right genus), and `21_2` `[1,14]`, `34_3` `[1,34]` (right genus,
wrong point counts: factors in mixed Hauptmodul coordinates).

**Regenerating.** The default recipe produces them since the fibre-product stage runs last in
`AllEquationsAboveCovers`; the pipeline's own run will write them in the same format, and from
2026-10-04 that stage never uses an equation the rebase or back-fill stage stored (those are in a
coordinate of their own). The writer here reads committed files, which do not record which entries
came from those stages; the trace-formula check is what protects it. ⚠ A fibre product of factors
from runs with DIFFERENT Hauptmodul normalisations can have the right genus and be the wrong curve
(`21_2 [1,14]`, `34_3 [1,34]` above); never accept one on genus alone.
## `21_2` and `6_5`: four quotients added 2026-10-03, both files regenerated with the default recipe

Both files predated the current cover propagation and were missing covers the pipeline now
produces; no flag is needed, the default recipe above fills them.

**`21_2`, the three genus-2 quotients `{1,3}`, `{1,6}`, `{1,14}`.** These are the quotients by the
three Atkin-Lehner involutions that act without fixed points, so Riemann-Hurwitz drops the genus
from 3 to 2 rather than to 1 or 0. They are now checked against **Guo-Yang's published curve**, not
against our own output: their pair is `y^2 = -(9u-1)(u+7)(u+3)`, `z^2 = -(u+3)`, `u = x^2`, with the
group acting by sign changes on `(x,y,z)`, and each quotient is the invariant subfield of one sign
pattern (the derivation is written out in `tests/X0_21_2.m`). Our models are those curves with `y`
scaled by `3/32`, `3/16`, `3/32`, so they agree over **Q**.

⚠ **The three were absent because the top curve here is stored as a non-hyperelliptic pair and the
old code only quotiented hyperelliptic curves** — but do not read that as the obstruction: a plain
rerun produces all three as hyperelliptic curves, so the limitation had already been lifted
elsewhere. The lesson is the usual one: rerun before diagnosing from the code.

**`6_5`, the quotient `{1,3,5,15}`.** The one genus-1 double cover of the star curve; its six
genus-0 siblings were already present. Gonzalez-Rotger publish quotients by a single `w_m`, so
nothing external covers this key; it rests on the independent checks in `ModelVerification.m`
(genus formula, L-polynomial divisibility over every nested pair, trace-formula point counts),
291 checks with no failures.

## Level M = 660: X_0^330(1), X_0^6(55), X_0^10(33), X_0^22(15) (2026-10-12, branch level-660-models)

The level-660 group (nine bases, 39 genus-0 quotients) was the largest single block behind the
Normaliz wall: its first polytope, (660, 231, 0), had run for seven days on lava in the old
coordinate order and died of its own time limit. With the coordinate order of PR #88 it took about
30 hours on lovelace at 8 threads, 12461 lattice points, and is committed here as
`polymake/polymake_solution_660_231_0`. The four even-D bases were built in the default recipe
(`genmodels.m`, BFCACHE=1, the composite-level tree with the exact m = 0 multipliers per support
class), with the solver script of PR #88 for any further rung and no weight-0 shift set. The five
odd-D bases, X_0^15(11), X_0^15(22), X_0^33(5), X_0^33(10) and X_0^55(3), need the two-cusp rung
(660, 68, 166), enumerated the same day (7591 points), and the odd-D fix of PR #91; they are running.

* `models_330_1.m`: 66 cover keys, 65 with an equation -- all 20 genus-0 quotients, covers up to
  genus 3, 6 stored as pairs over a conic -- and 1 empty. 8 minutes.
* `models_6_55.m`: 65 keys, 41 with an equation -- all 5 genus-0 quotients, covers up to genus 5,
  3 pairs over a conic -- and 24 empty. 22 minutes.
* `models_10_33.m`: 65 keys, 41 with an equation -- all 5 genus-0 quotients, covers up to genus 7,
  3 pairs over a conic -- and 24 empty. 22 minutes.
* `models_22_15.m`: 65 keys, 28 with an equation -- both genus-0 quotients, covers up to genus 6 --
  and 37 empty. 25 minutes.

Outside check: every stored hyperelliptic curve of genus >= 1 has the point counts of the
Eichler-Selberg trace formula over F_p and F_{p^2} at the six smallest primes not dividing 330,
p = 7, 13, 17, 19, 23, 29 (campaign `composite/checkcomp.m`): 234 of 234, 198 of 198, 198 of 198
and 156 of 156. No published equation exists for any of these bases.
## Level M = 1020: X_0^510(1), X_0^6(85), X_0^10(51) (2026-10-11, branch level-1020-models)

The first polytope every base of the level-1020 group asks for, (1020, 360, 0), is committed here
as `polymake/polymake_solution_1020_360_0` (11973 lattice points), enumerated on lovelace in
68 minutes with the coordinate order of PR #88 where the old order had not finished after five
days. The three even-D bases were built in the default recipe (`genmodels.m`, BFCACHE=1, the
composite-level tree with the exact m = 0 multipliers per support class), with the solver script of
PR #88 for any further rung and no weight-0 shift set. The three odd-D bases, X_0^15(17),
X_0^15(34) and X_0^51(5), wait on the two-cusp rungs.

* `models_510_1.m`: 65 cover keys, 55 with an equation -- all 11 genus-0 quotients, covers up to
  genus 5 -- and 10 empty. 33 minutes.
* `models_6_85.m`: 65 keys, 28 with an equation -- both genus-0 quotients, covers up to genus 5 --
  and 37 empty. 61 minutes.
* `models_10_51.m`: 50 keys, 22 with an equation -- the one genus-0 quotient, covers up to genus 6 --
  and 28 empty. 59 minutes.

Outside check: every stored curve of genus >= 1 agrees with the trace formula over F_p and F_{p^2}
at the six smallest primes not dividing 510, p = 7, 11, 13, 19, 23, 29: 264 of 264, 156 of 156 and
126 of 126. No published equation exists for any of these bases.
## Level M = 780: X_0^390(1), X_0^6(65), X_0^10(39), X_0^26(15) (2026-10-11, branch level-780-models)

The first polytope every base of the level-780 group asks for, (780, 274, 0), is committed here as
`polymake/polymake_solution_780_274_0` (12729 lattice points). It was enumerated on lovelace in
74 minutes with the coordinate order of PR #88 (divisors farthest from sqrt(M) first, Normaliz
keeping the order); the run in the old order had not finished after five days. The four even-D
bases were then built in the default recipe (`genmodels.m`, BFCACHE=1, the composite-level
tree with the exact m = 0 multipliers per support class), each in 16 to 35 minutes, with the solver
script of PR #88 for any further rung. No weight-0 shift set was needed. The two odd-D bases,
X_0^15(13) and X_0^39(5), wait on the two-cusp rungs, as at level 924.

* `models_390_1.m`: 66 cover keys, 55 with an equation -- all 11 genus-0 quotients, covers up to
  genus 5 -- and 11 empty. 16 minutes.
* `models_6_65.m`: 65 keys, 41 with an equation -- all 5 genus-0 quotients, covers up to genus 5 --
  and 24 empty. 33 minutes.
* `models_10_39.m`: 65 keys, 33 with an equation -- all 3 genus-0 quotients, covers up to genus 4,
  5 of them stored as pairs over a conic -- and 32 empty. 33 minutes.
* `models_26_15.m`: 50 keys, 22 with an equation -- the one genus-0 quotient, covers up to genus 6 --
  and 28 empty. 35 minutes.

Outside check: every stored hyperelliptic curve of genus >= 1 has the point counts of the
Eichler-Selberg trace formula over F_p and F_{p^2} at the six smallest primes not dividing 390,
p = 7, 11, 17, 19, 23, 29 (campaign `composite/checkcomp.m`): 264 of 264, 216 of 216, 150 of 150
and 126 of 126. No published equation exists for any of these bases.
## Level M = 924: X_0^462(1) (2026-10-06, branch level-924-models)

The first base of the level-924 group (eight bases, 21 genus-0 quotients). Its Borcherds search asks
first for the polytope (924, 325, 0), enumerated once on lava (Normaliz, 8 threads, 4 h 19 min, 8674
lattice points) and committed here as `polymake/polymake_solution_924_325_0`, together with the
weight-0 shift set `polymake/tshift_w0_924.txt` (283 points, minutes) that the t-shift fallback uses
for the taller m = 0 rungs; the cache is keyed by the polytope alone, so every level-924 base shares
them. `models_462_1.m`: default recipe (`genmodels.m`, BFCACHE=1, the composite-level tree at commit
e85261a of branch exact-m0, PR #79, which carries the m = 0 term of #66, #74 and #79) on lovelace,
15 minutes; 66 cover keys, 55 with a model (49 as y^2 = f(x), 6 as double covers of a conic) -- all 11
genus-0 quotients of the base, covers up to genus 5 -- and 11 empty. Outside
check: every stored curve of genus >= 1 has the point counts of the Eichler-Selberg trace formula
over F_p and F_{p^2} at the six smallest primes not dividing 462, p = 13, 17, 19, 23, 29, 31 (228 of
228; campaign `composite/checkcomp.m`; the six double covers of a conic were not part of this check). No published equation exists for this base.

`models_6_77.m` (2026-10-06): the first composite-level base of the group (N = 77), same recipe on
lovelace, 26 minutes, the m = 0 multipliers from the algebraic routine per support class {7}, {11},
{7, 11}. 65 cover keys, 37 with a model (34 as y^2 = f(x), 3 as double covers of a conic) -- all 4 genus-0
quotients, covers up to genus 7 -- and 28 empty. Outside check: every stored curve of genus >= 1 agrees with the
trace formula over F_p and F_{p^2} at the six smallest primes not dividing 462 (180 of 180; the three
double covers of a conic were not part of this check). At
p = 11, a prime of the level, the genus-1 quotient by {1, 2, 7, 14, 33, 66, 231, 462} reduces to a
smooth curve whose counts differ from the trace formula; the trace formula does not describe the
reduction at a prime of bad reduction of the Shimura curve, so the check now skips such primes. No
published equation exists for this base.

`models_22_21.m` (2026-10-09): same recipe on lovelace, 34 minutes, the m = 0 multipliers per
support class {3}, {7}, {3, 7}. 50 cover keys, 22 with an equation -- the one genus-0 quotient
{1, 3, 14, 22, 42, 66, 77, 231} (the star curve is the Hauptmodul base and is stored nowhere),
covers up to genus 6 -- and 28 empty, none of which the Hauptmodul rebase sweep fills. Outside
check: every stored curve of genus >= 1 agrees with the trace formula
over F_p and F_{p^2} at the six smallest primes not dividing 462 (126 of 126). No published equation
exists for this base.

The other four bases with N > 1 of the group, X_0^21(11), X_0^21(22), X_0^33(7) and X_0^33(14), are
NOT here. After the {oo}-weakly holomorphic ring (the shared rung above), their Borcherds searches
go on to ask for the ring of forms with poles at both cusps 0 and oo, which the three bases above
never request. That search climbs the weight-236 ladder (924, n, 236) in steps of 236 in n, finds
rank 0 at every rung while the requested dimension grows by 236, and dies at the cap n = 10000 after
3.5 days. The climb depends on the level alone: the logs of the four runs are byte-for-byte identical
up to the point each was stopped.

**Every empty key of `22_21` and `6_77` filled as a fibre product (2026-10-11, branch
level-924-fibre-products).** The 28 empty keys of each file are the covers with no genus-0
Atkin-Lehner quotient to carry a hyperelliptic equation; each is the fibre product over the star
line of its index-2 double covers already in the file, built by the writer of the fibre-product
stage (`tools/writefp.m`, as the 203 entries above) and stored in the same `<genus, "CRV", [...]>`
shape: `22_21` genera 3 to 7, `6_77` genera 3 to 9. The writer accepted each one on the genus formula
and the trace formula over F_p and F_{p^2} at the first three good primes. Then, independently of
the writer, every entry was re-read from its stored strings alone and compared with the trace
formula at the 4th and 5th good primes, which the writer never used: 112 of 112 (prime, cover)
counts agree. `462_1` keeps its 11 empty keys: none of them has committed double covers to multiply.
As for the 203 above, at level N > 1 these rest on the m = 0 term of Schofer's formula through their
factors (PR #66), so they are checked by point counts, not proved.
The five odd-D bases of the group, X_0^21(11), X_0^21(22), X_0^33(7), X_0^33(14) and X_0^77(3), are
not here. For odd D the search also needs the ring of weight-1/2 forms with poles at both cusps,
which asks for the two-cusp polytope (924, 89, 236). That enumeration did not finish within the
two-hour Normaliz limit the runs used, and a timed-out solve then returned "no points" rather than
an error (PR #81 makes it an error), so the four runs that were started climbed the ladder
(924, n, 236) at rank 0 for 3.5 days until the cap n = 10000. The two-cusp rung has since been
enumerated (3194 points, with the coordinate order of PR #88); the bases with N > 2 also need the
divisor-constraint fix of PR #91, and the odd-D bases follow in a later PR.
## The odd-D bases behind the wall, first base: X_0^33(10) (2026-10-13, branch odd-d-models-2026-10-13)

The odd-D bases of the level groups 660, 780, 924 and 1020 need, besides the shared rungs with poles
at oo only, the ring of weight-1/2 forms with poles at both cusps, and that search was stopped twice:
once by the two-cusp polytopes, which Normaliz did not enumerate within the runs' time limits
(PRs #80 and #88), and once by the pool of forms at the cusp 0, which grows with the product of the
cusp-0 order and the level and took the kernel computation of PR #63 to a 172 GB segmentation fault at
X_0^33(7) (PR #95 replaced it by a forward substitution against the triangular pool). Two further
defects affected only odd D with N > 2: the divisor of a Borcherds form was sampled at the wrong
discriminants (PR #91), and the expansion at the cusp 0 was read in the wrong variable (PR #63).
With the three fixes merged and the exact m = 0 multiplier of PR #79 (not merged) applied on top, the
fifteen bases were started together on lovelace on 2026-10-13 from one tree (`genmodels.m`,
BFCACHE=1, 8 threads, 10-day cap); reproducing these files from main therefore needs PR #79, and
PR #105 for the fibre products. This section records them as they finish.

**`models_33_10.m`** (31 minutes). 50 cover keys, all with an equation: 22 from the Borcherds run
-- the one genus-0 quotient, by {1, 2, 15, 30, 33, 55, 66, 110}, is y^2 = (88209/16) x (1 - x), a
P^1 with a rational point; then 5 of genus 1, 7 of genus 2, 3 of genus 3, 3 of genus 4, 2 of genus 5
and one of genus 8 -- and 28 fibre products of the file's double covers (the paragraph below). Its
searches asked for the polytopes (660, 231, 0) (already committed with level 660), (660, 401, 0) and
the two-cusp polytope (660, 69, 166); the last two are committed here. `polymake_solution_660_69_166`
is a Normaliz enumeration on lovelace, 8170 lattice points (as is the committed (660, 231, 0),
12461 points; in each file the first line is the Magma sequence header and the last its closing
bracket). `polymake/polymake_solution_660_401_0` is NOT the lattice points of its polytope: the
search was served by the shifted-products fallback of `nmzsolve.py`, which built a spanning set of
855 points from the rung (660, 231, 0) and the 441 weight-0 shifts (the run's log line: "base rung
231 (12461 pts) + 441 shift(s) x1 -> 855 pts for n=401"); the file carries the first-line mark
`// derived: ...` that the fallback writes, so the lower-rung read and the products path do not
take it for a complete enumeration. With these files the default recipe reproduces the base
without Normaliz. Outside check: every stored curve of genus >= 1 has the point counts of the
Eichler-Selberg trace formula over F_p and F_{p^2} at the six smallest primes not dividing 330,
p = 7, 13, 17, 19, 23, 29 (campaign `composite/checkcomp.m`): 126 of 126. Six primes do not by
themselves determine the quadratic twist: several quotients of genus 2 to 4 here have trace zero at
every one of these primes where the character of a twist by -33, -3, -30, 231, -770 or 2310 is -1,
so the twisted curve has the same six counts. Primes up to 113 rule every such twist out (Sachi's
independent check at the 25 primes from 13 to 113: all counts agree). `tools/check_models_trace.m`
(this PR) compares all 49 curves of positive genus, hyperelliptic entries and fibre products alike,
with the trace formula over F_p and F_{p^2}: at the 13 primes from 13 to 61, 637 of 637 prime
counts agree (each prime counting both fields), none skipped. The run with the same script at the primes to 113 was killed for memory on the
Mac (two Magma processes at once) and is to be repeated on lovelace. No published equation
exists for this base. As at every base with N > 1, the equations rest on the m = 0 term of Schofer's
formula (PRs #66, #74, #79, not yet merged), and nothing fixes which quadratic twist of the conic the
genus-0 quotient is beyond the covers above it agreeing with one another.

**`models_33_14.m`** (12 h 07 min, 16 threads; the deep rung is the pool of 7624 forms at the cusp 0
at order m = -33, whose good forms took 950 s of forward substitution -- the step that had been the
172 GB kernel). 50 cover keys, all with an equation: 15 from the Borcherds run -- 5 of genus 1, 2 of
genus 2, 3 of genus 3, 2 of genus 4, 3 of genus 5 -- and 35 fibre products of them (the paragraph
below); this base has no genus-0 quotient other than the star curve itself, so every equation here
is a cover of the Hauptmodul line. Its searches asked for (924, 325, 0) (committed with level 924)
and the two-cusp polytope (924, 89, 236), committed here (`polymake/polymake_solution_924_89_236`,
3194 lattice points, the Normaliz enumeration of 2026-10-12 on lovelace that PR #80 had recorded as
missing). Outside check: every stored curve has the point counts of the trace formula over F_p and
F_{p^2} at p = 13, 17, 19, 23, 29, 31 (campaign `composite/checkcomp.m`): 90 of 90; and at the 25
primes from 13 to 113, which separate every quadratic twist by a d supported on 2, 3, 7, 11 and -1,
`tools/check_models_trace.m` compares all 50 curves of positive genus with the trace formula over
F_p and F_{p^2}: at the 13 primes from 13 to 61, 650 of 650 prime counts agree, none skipped
(Sachi's independent check at the 25 primes from 13 to 113 agrees throughout; the run to 113 with
this script was killed for memory on the Mac and is to be repeated on lovelace). No published
equation exists for this base; the same remarks on the m = 0 term apply.

**Every empty key of both files filled as a fibre product (2026-10-14).** The 28 empty keys of
`33_10` and the 35 of `33_14` are the covers with no genus-0 Atkin-Lehner quotient to carry an
equation; each is the fibre product over the star line of its index-2 double covers already in the
file, built by the writer of the fibre-product stage (`tools/writefp.m`) with the coordinate check of
PR #105: for every generating set the three double covers above each pair of factors were tested,
and at all 63 the product of the three is a square, so the factors are consistent with one
coordinate on the star line (up to a quadratic twist shared by two of them, which the test cannot
see). Genera 3 to 8 at `33_10`
(2, 2, 5, 4, 9, 6 of genus 3 through 8) and 3 to 11 at `33_14` (2, 4, 1, 0, 3, 4, 8, 8, 5 of genus 3
through 11, so `33_14`'s top curve of genus 11 is among them). The writer accepted each on the genus
formula and the trace formula over F_p and F_{p^2} at the first three good primes; then,
independently of the writer, every entry was re-read from its stored strings and compared with the
trace formula at the 4th and 5th good primes: 126 of 126 counts agree. Stored in the
`<genus, "CRV", [...]>` shape, with the note at the top of each file.

## The reachable odd-D backlog, first two bases (2026-10-05, branch backlog-odd-2026-10)

`models_33_2.m` (15 cover keys: 2 of genus 0, 6 of genus 1, 4 of genus 2, 3 of genus 3) and
`models_93_2.m` (12 keys: 1 of genus 0, 2 of genus 1, 1 of genus 2, 2 of genus 4, 1 of genus 8),
default recipe (`genmodels.m`, BFCACHE=1) on lava from branch composite-level, 3 min and 46 min.
Outside check: every stored curve of genus >= 1 has the point counts of the Eichler-Selberg trace
formula over F_p and F_{p^2} at p = 11, 13, 17, 19, 23, 29 (50 of 50 and 36 of 36; campaign
`composite/checkcomp.m`). No published equation exists for either base. `models_55_2.m` (15 keys: 1 of genus 0, 3 of genus 1, 2 of genus 2, 3 of
genus 3, 1 of genus 6; 68 min) and `models_87_2.m` (14 keys, genus up to 8; 2 h): their first runs
stopped in the exact m = 0 multiplier's class-constancy check at d = -15, a roundoff miss at Prec 80
(3e-15 against a floor of 1e-15 at M = 220, and 2e-5 at M = 348); with the precision retry of
composite-level 9201ae0 both passed at Prec 160 and the trace formula agrees at every curve (45 of 45
each). `models_69_2.m` (14 keys: 3 of genus 1, 4 of genus 2; 7 keys stay empty after the Hauptmodul
sweep, none of them genus 0) ran the same way on lava and agrees with the trace formula at every
curve (35 of 35). This base has no genus-0 quotient, so its equations are new covers only.
`models_95_2.m` (14 keys: 1 of genus 0, 2 of genus 1, then genus 2, 3, 4, 5, 5, 7, 10; 4 keys empty)
ran on lava the same way (12 h; its m = 0 multipliers missed the Prec-80 floor once and passed at
Prec 160) and agrees with the trace formula at every curve of genus >= 1 (45 of 45).
`models_111_2.m` (12 keys: 1 of genus 0, 2 of genus 1, then genus 2, 5, 5, 10; 5 keys empty): its first
run (composite-level code) lost all precision in the numerical m = 0 multiplier at M = 444 (a
class-constancy deviation of 3e61) and sat in the Prec-160 retry for 13 hours; restarted from the
cached Borcherds forms with the exact algebraic multiplier (branch exact-m0), it finished in 25 minutes.
Trace formula 36 of 36.

## Three stale files rebuilt: X_0^22(13), X_0^34(5), X_0^38(7) (2026-10-06)

The old files held only the quotients by a single involution (most of them empty) and the genus-1 and
genus-2 covers, 5 to 10 keys each, with no provenance; their one genus-0 quotient each -- by the
order-4 groups {1, 11, 26, 286}, {1, 5, 34, 170}, {1, 7, 38, 266} -- was absent, which a census had
read as a conic without a rational point. Rerun with the current pipeline (default recipe, branch
exact-m0, 4 to 8 minutes each on a laptop) they hold 14, 15 and 13 keys, with the genus-0 quotient
present and rational in every case: y^2 = (121/16) x (x - 1), y^2 = (2x - 3)(2x + 1) and y^2 = 5 - 8x,
each with an evident rational point, so all three curves are P^1 over Q and nothing here was a conic.
Covers up to genus 5 (22_13), 4 (34_5) and 4 (38_7). Trace formula at every curve of genus >= 1:
36 of 36, 45 of 45, 35 of 35 (campaign composite/checkcomp.m). No published equation covers these
keys. The lesson is the standing one: a missing quotient is a stale file until a rerun says otherwise.

Three keys in these files are hyperelliptic quotients that are still empty and have not yet been
computed, as opposed to being out of the method's reach: `93_2` by w_186 (`{1,186}`), `111_2` by w_222
(`{1,222}`) and `38_7` by w_7 (`{1,7}`). PR #73 stores models for the quotients by w_2 of `22_13` and
`38_7`, as pairs of equations over a conic, which are empty here; the models the two branches share
agree up to isomorphism.
## `models_210_1.m` and the level-420 polytope cache (2026-10-05, branch level-420-models)

`X_0^210(1)` was behind the Normaliz wall (M = 420 has 24 divisors) until its first polytope rung
was enumerated in August for the theta campaign: `polymake/polymake_solution_420_145_0` (19864
lattice points, 13.4 wall-hours of Normaliz, banked on the campaign branch as
`vvdata/weyl-campaign/nmz_420_145_0.txt.gz`). The t-shift fallback of `nmzsolve.py` answers the
second rung from it (`polymake_solution_420_250_0`, 212 points, a spanning set of the rung's forms and
not the lattice points of the polytope; its first line says so, as a Magma comment, so that a later
enumeration is not mistaken for it); both files are committed here so the default recipe reproduces
the base without Normaliz. The run used the composite-level tree on lovelace at ca87e19 (branch
composite-level, PR #74; at N = 1 no m = 0 term is involved). With them the whole run took 7 min 37 s on
lovelace: 66 cover keys, 40 of genus 0, 31 of genus 1, 10 of genus 2, 4 of genus 3 and 6 fibre
products of genus 3. Outside check: 234 of 234 point counts of the Eichler-Selberg trace formula at
p = 11, 13, 17, 19, 23, 29 (campaign `composite/checkcomp.m`). The same cache serves every base with
M = 420: the composite-level bases `6_35`, `10_21`, `14_15` (branch composite-level) ran from it; the
bases that also need a cusp-0 rung, `15_7`, `21_5`, `15_14`, `21_10`, `35_6` (all `(420, 42, 104)`)
and `35_3` (`(420, 146, 104)`), wait on those enumerations.
## 36 quotient models from a rerun of the gap bases (2026-10-03, lovelace, code at e4d8b4a, a commit on the branch of PR #63 and not on `main`)

The census of 2026-10-03 found 100 subhyperelliptic quotients with no model across 31 bases, almost
all of them because the committed file predates a code fix rather than because the method fails
(`21_2` and `6_5` closed by a plain rerun, PR #70). The 27 reachable bases were rerun on lovelace
with the default recipe (`genmodels.m`, `NORMALIZ_BIN=/usr/bin/normaliz`, 10 h cap, all 24 that
finished exited 0; `143_1`, `145_1`, `203_1` run separately). Against this branch the reruns add
36 entries in 8 files, of genus 0 to 7 (3, 8, 9, 6, 5, 3, 1, 1 entries of genus 0 through 7):

    14_19  [1,133] [1,7,19,133] [1,19] [1,14,38,133]
    22_13  [1,286] [1,26] [1,22,26,143] [1,13,22,286] [1,11,26,286] [1,11,13,143] [1,11]
    22_7   [1,22]            26_5  [1,10]
    34_5   nine: [1,5,34,170] [1,10,34,85] [1,2,5,10] [1,5] [1,34] [1,5,17,85] [1,170] [1,2,85,170] [1,10,17,170]
    34_7   [1,14,34,119] [1,7,34,238]
    38_7   [1,14,19,266] [1,2,7,14] [1,266] [1,7,38,266] [1,38] [1,14,38,133]
    58_5   [1,10] [1,10,58,145] [1,2,5,10] [1,29] [1,2,29,58] [1,5,29,145]

Where a key already held a model of the form y^2 = f(x) the rerun's version (often a different
normalisation of the same curve) was discarded. Where this branch held a fibre-product `CRV` pair
and the rerun found a model y^2 = f(x) for the same key, the hyperelliptic model replaces the pair
(Eran's decision, 2026-10-13: the hyperelliptic form is the more useful one when it exists; the
two were produced independently and both agree with the trace formula): 13 keys, `10_17 [1,85]`,
`10_29 [1,58]`, `10_37 [1,2] [1,185]`, `14_13 [1,91]`, `14_19 [1,14]`, `22_7 [1,11]`, `26_5 [1,65]`,
`34_3 [1,3] [1,6]`, `34_7 [1,14] [1,238]`, `6_41 [1,123]`, each with the genus of the pair it replaces
(the rerun output is in `gapscan/out` on lovelace). The two
odd-$D$ reruns `115_1` and `161_1` reproduce the entries of PR #69 exactly. Vouching: the same
checks as every pipeline model (`tests/ModelChecks`, genus, and the published tables where they
reach: `GonzalezRotger`, `GuoYangTable1`).

## Twenty quotients from the heights of their branch points (branch models-1326-210)

Every quotient in this section has genus 3 and lies behind the divisor wall of the Borcherds
search, so none of the recipes above reaches it. Each has an elliptic Atkin–Lehner quotient `E =
X/<W, w>` and is a double cover of it branched at four CM points. The heights of those points on
`E` are given by the explicit Gross–Zagier formula of Cai–Shu–Tian (Algebra & Number Theory 8
(2014), Theorem 1.5, second special case), implemented in `HeegnerHeights.m`; the points of `E` of
those heights form a finite set, and every configuration is tested exactly against the trace
formula by `ExactDoubleCover.m`. Nothing in the construction is numerical except the L-values,
which direct the search; what certifies a result is the agreement of the point counts. The files
hold the genus-3 quotient only.

**What the formula gives, and what is enumerated.** For a branch place of degree 2 over a
quadratic field `F`, the character attached to `F` gives the height of `P - Pbar` exactly, so `P -
Pbar` is a known multiple of the generator of `E(F)^-`; that multiple being an integer is also the
test that identifies `E` inside its isogeny class, since on the other curves it is not. The
rational part is the trivial character's component, which gives the sum of the points over a place
of degree `e` as `e·xi_E ± b·G`, with `xi_E` the image of the Hodge class of `X_0^D(N)`, common to
every place, and `b` a rational number. So the enumeration runs over the signs and over the
multiple of the generator in `xi_E`, which matters only modulo 1 because translating every branch
point by a rational point of `E` gives an isomorphic cover, and which has denominator dividing
twice the denominators of the `b`'s. That makes the enumeration finite and complete. When every
branch place has degree 1 there is no minus part, hence no lattice test, and the integrality of
the rational parts takes its place.

**Quotients whose elliptic factor is an oldform.** At five of these the quaternionic form is old
at 2, so `E` has conductor `DN/2` and the routine that reads the quotient off the Atkin–Lehner
signs does not find it; `E` was identified from the trace formula, by matching the point counts of
the genus-one quotient against the elliptic curves of conductor dividing `DN`. The formula then
needs the degree of the map it is applied to. Write `p = 2`, let `f` parametrise `E` from the
curve on which its form is new, and let `pi_1`, `pi_2` be the two degeneracy maps, `pi_2 = pi_1
w_p`. The quotient map factors through `w_p`, so it is `g = f pi_1 + f pi_2`, of degree `2(p + 1 +
a_p)·deg f`: each `f pi_i` has degree `(p + 1) deg f` and each cross term is `f T_p f^dual`, of
degree `a_p deg f`. Where `p` ramifies in the field of the CM point, `w_p` acts on it as the
Frobenius of the prime above `p`, so `pi_2(P)` is that conjugate of `pi_1(P)` and the
`chi`-component of `g(P)` is `(1 + chi(p))` times the component of `f(pi_1 P)`. With the degree
taken to be `deg g`, every component of the formula is therefore multiplied by

    (1 + chi(p))^2 / (2(p + 1 + a_p)),

which is `2/(p + 1 + a_p)` for the characters with `chi(p) = +1` and zero for the others. That is
2 at 141a1, 357d1 and 345f1, where `a_2 = -2`, and 1 at 285a1, where `a_2 = -1`. The vanishing at
`chi(p) = -1` is what explains the rational branch places whose nontrivial characters would
otherwise have had nonzero components. The argument is for `p` ramified in the field of the branch
point; the same expression is used at a place where it is not, which is not derived, so there the
height directs the search and only the counts certify the answer. This derivation is Sachi
Hashimoto's.

**Reproducing.** `tools/heegner_model.m` takes the base, the involution, and for an oldform
quotient the curve and the prime at which the form is old:

    magma -b D:=1326 N:=1    W:=1,17,78,1326 U:=6 tools/heegner_model.m < /dev/null
    magma -b D:=210  N:=13   W:=1,2,21,42,65,130,1365,2730 U:=5 tools/heegner_model.m < /dev/null
    magma -b D:=39   N:=11   W:=1,11,39,429 U:=3 tools/heegner_model.m < /dev/null
    magma -b D:=55   N:=7    W:=1,35,55,77 U:=5 tools/heegner_model.m < /dev/null
    magma -b D:=6    N:=301  W:=1,3,14,42,43,129,602,1806 U:=903 tools/heegner_model.m < /dev/null
    magma -b D:=6    N:=455  W:=1,6,7,10,13,15,42,70,78,91,105,130,195,546,910,1365 U:=30 tools/heegner_model.m < /dev/null
    magma -b D:=6    N:=187  W:=1,3,11,17,33,51,187,561 U:=374 tools/heegner_model.m < /dev/null
    magma -b D:=22   N:=51   W:=1,6,17,22,33,102,374,561 U:=187 tools/heegner_model.m < /dev/null
    magma -b D:=6    N:=247  W:=1,2,19,38,39,78,741,1482 U:=3 tools/heegner_model.m < /dev/null
    magma -b D:=6    N:=253  W:=1,2,23,33,46,66,759,1518 U:=3 tools/heegner_model.m < /dev/null
    magma -b D:=6    N:=185  W:=1,2,15,30,111,185,222,370 U:=3 tools/heegner_model.m < /dev/null
    magma -b D:=6    N:=385  W:=1,5,6,11,14,21,30,55,66,70,105,154,231,330,770,1155 U:=2 tools/heegner_model.m < /dev/null
    magma -b D:=714  N:=5    W:=1,3,10,14,17,30,35,42,51,105,170,238,510,595,714,1785 U:=2 tools/heegner_model.m < /dev/null
    magma -b D:=15   N:=38   W:=1,2,15,19,30,38,285,570 U:=2 tools/heegner_model.m < /dev/null
    magma -b D:=141  N:=2    W:=1,3,94,282 U:=2 curve:=141a1 old:=2 tools/heegner_model.m < /dev/null
    magma -b D:=21   N:=34   W:=1,3,14,17,42,51,238,714 U:=2 curve:=357d1 old:=2 tools/heegner_model.m < /dev/null
    magma -b D:=57   N:=10   W:=1,5,6,19,30,95,114,570 U:=2 curve:=285a1 old:=2 tools/heegner_model.m < /dev/null
    magma -b D:=15   N:=46   W:=1,5,6,23,30,115,138,690 U:=2 curve:=345f1 old:=2 tools/heegner_model.m < /dev/null
    magma -b D:=690  N:=7    W:=1,6,10,14,15,21,23,35,138,210,230,322,345,483,805,4830 U:=2 tools/heegner_model.m < /dev/null

**The models.** `E` is the elliptic quotient, and the branch data gives the CM discriminant of
each branch place with the field it is defined over.

* **X_0^1326(1)/{1,17,78,1326}: `y^2 = 272x^8 - 72x^6 - 55x^4 - 18x^2 + 17`** (stored reversed),
  over 1326b2, rank 1, branched at discriminants −24 and −408 over `Q(sqrt 2)`. The isogeny class
  1326b is given by the Atkin–Lehner signs and the heights pick b2 within it. Obtained once; the
  surviving configuration has both sums equal to the generator.
* **X_0^210(13)/{1,2,21,42,65,130,1365,2730}: `y^2 = 4x^8 - 20x^7 + 97x^6 - 20x^5 - 138x^4 + 20x^3
  + 97x^2 + 20x + 4`**, obtained twice, over 2730s1 (branch points −312, −120 over `Q(sqrt 2)`)
  and over 2730u1 (−420, −1092 over `Q(i)`), both of rank 1; the two curves are Q-isomorphic.
* **X_0^39(11)/{1,11,39,429}: `y^2 = -891x^8 + 1044x^6 + 214x^4 - 124x^2 + 13`**, over 429b1 (the
  lattice test excludes the other five curves of 429b). Branch points −132 over `Q(sqrt -11)` and
  −52 over `Q(sqrt 13)`, with the genus characters (−11, 12) and (−4, 13).
* **X_0^55(7)/{1,35,55,77}: `y^2 = -35x^8 - 2028x^6 + 198x^4 + 68x^2 + 5`**, over 385a1. Branch
  points −1540 over `Q(sqrt -35)` and −20 over `Q(sqrt 5)`, characters (−35, 44) and (−4, 5).
* **X_0^6(301)/{1,3,14,42,43,129,602,1806}: `y^2 = 129x^8 - 2276x^6 + 3334x^4 + 1180x^2 - 63`**,
  over 1806l1, the only curve of its class; branch points −1204 over `Q(sqrt -7)` and −1032 over
  `Q(sqrt 129)`, minus parts four times the generators.
* **X_0^6(455)/W: `y^2 = -416x^8 - 515x^6 + 582x^4 - 251x^2 + 40`**, over 2730y1; branch points
  −10920 over `Q(sqrt -26)` and −1560 over `Q(sqrt 10)`. 2730y2 is excluded by the lattice test,
  and 2730y3, 2730y4 pass it but no configuration over them is consistent with the trace formula,
  which is the expected behaviour: the lattice test is necessary, not sufficient.
* **X_0^6(247)/W: `y^2 = 81x^8 - 4788x^6 + 1654x^4 + 2540x^2 + 513`**, over 1482i1; the rational
  CM points of discriminants −52 and −3 and a place of degree 2 of discriminant −456 over `Q(sqrt
  57)`.
* **X_0^6(253)/W: `y^2 = 9x^8 + 1388x^6 + 3526x^4 - 3060x^2 - 1863`**, over 1518p1; discriminant
  −552 over `Q(sqrt -23)` and the rational points −88 and −1012. Here the multiple of the
  generator in `xi_E` has denominator 3, which is why an earlier version of the search that
  assumed it integral found nothing.
* **X_0^6(185)/W: `y^2 = 2025x^8 + 4716x^6 - 122x^4 - 5620x^2 - 999`**, over 1110m1; the rational
  point −555, a place of degree 2 of discriminant −4440 over `Q(sqrt -111)`, and the rational
  point −40. The numerical fit of Nelson's method here had recognised only one of its three
  invariants, so this is the first model of that quotient.
* **X_0^6(385)/W: `y^2 = 9x^8 - 580x^6 - 9482x^4 + 60092x^2 + 3465`**, over 2310q1; the rational
  points −40 and −1540 and a place of degree 2 of discriminant −9240. See the caveat below.
* **X_0^714(5)/W: `y^2 = 9720x^8 + 2061x^6 + 598x^4 - 155x^2 + 16`**, over 3570x2; a place of
  degree 2 of discriminant −14280 and the rational points −840 and −340. Two configurations
  survive, with twists −43 and −2171, and they are Q-isomorphic. The search needs `Effort := 2` in
  the Mordell–Weil computations. See the caveat below.
* **X_0^15(38)/{1,2,15,19,30,38,285,570}: `y^2 = 81x^8 + 180x^6 + 2062x^4 - 1660x^2 + 361`**, over
  570e4. Every branch place is rational, of discriminants −228, −12, −760 and −40, so there is no
  lattice test; the integrality of the rational parts takes its place and excludes the other three
  curves of 570e (on 570e1 and 570e3 the parts come out 1, 3/2, 1/2, which no common shift makes
  integral, and on 570e2 they are not rational). The branch point of discriminant −12 belongs to
  an order of conductor 2, which the formula as implemented does not reach, so its position was
  searched over the integers of absolute value at most the window rather than predicted, and came
  out at 0, that is, the point is `xi_E` up to torsion. **Uniqueness there rests on that window.**
* **X_0^141(2)/{1,3,94,282}: `y^2 = 81x^8 + 252x^6 - 26x^4 - 260x^2 - 47`**, over 141a1, form old
  at 2; the rational points −4 and −24 and a place of degree 2 of discriminant −564 over `Q(sqrt
  -47)`.
* **X_0^21(34)/{1,3,14,17,42,51,238,714}: `y^2 = 81x^8 - 612x^6 + 262x^4 + 220x^2 + 49`**, over
  357d1, form old at 2.
* **X_0^57(10)/{1,5,6,19,30,95,114,570}: `y^2 = 729x^8 + 500x^6 - 2x^4 + 68x^2 - 15`**, over
  285a1, form old at 2, where the multiplier is 1 rather than 2.
* **X_0^15(46)/{1,5,6,23,30,115,138,690}: `y^2 = 81x^8 - 900x^6 + 358x^4 + 380x^2 + 81`**, over
  345f1, form old at 2; one of its places sits at ten times the generator, outside the search
  window the earlier version used.
* **X_0^690(7)/W: `y^2 = 9x^8 - 130x^6 - 767x^4 - 520x^2 + 144`**, over 4830bj1, branched at four
  rational points, of discriminants −3, −280, −420 and −115. This base was out of reach until the
  search for an optimal embedding of a CM order was corrected: it reported that the order of
  discriminant −280 does not embed in the Eichler order of level 7 in `B_690`, although it does, so
  no branch data could be computed at all.
* The remaining two are **hyperelliptic with no model `y^2 = f(x)` over Q**: their canonical image
  is the conic ramified at 2 and at the real place, which has no rational point to project from.
  They are stored as double covers of that conic, in the shape of the other paired entries in this
  directory, with coordinates `s`, `z` of weight 1 and `y` of weight 2, and `x` the third
  coordinate of the conic scaled so that the equation is integral: **X_0^6(187)/W** is `X1^2 +
  2X2^2 + 22X3^2 = 0` with `y^2 = 4449s^4 + 84056s^3z + 411472s^2z^2 + 195328sz^3 + 797696z^4`,
  over 1122i1 (1122i2 excluded by the lattice test), branched at −4488 and −264 over `Q(sqrt -2)`;
  two configurations survive, with twists −14 and 721, and they give the same curve, the
  isomorphism being exhibited rather than inferred from the counts. **X_0^22(51)/W** is `X1^2 +
  2X2^2 + 59X3^2 = 0` with `y^2 = 9242176s^4 - 2871136s^3z + 4331508s^2z^2 + 124292sz^3 +
  4787449z^4`, over 1122f1, branched at the same two discriminants, one surviving configuration,
  twist −7. The 59 is an artefact of the coordinates and not bad reduction: that conic has points
  over every `Q_p` with `p` odd.

**X_0^6(323)/{1,6,19,34,51,114,646,969}, where the branch locus is one place of degree 4.** The four
points are a single Galois orbit, so the characters that would pin one of them have order 4, and
their L-functions do not factor into quadratic twists: the formula as implemented does not reach
them. What it does give is enough to leave a small search. The trivial character puts the sum of the
four points at `4·xi_E ± 2G`; the character attached to the quadratic subfield `Q(sqrt -2)` of the
place's field has a vanishing component, because the twist of 1938h1 by −8 has analytic rank 2, so
`P + sigma^2 P` and its conjugate differ by torsion alone, with `sigma` a generator of the cyclic
group of order four acting on the orbit. Writing `F` for the quartic field (discriminant 62016) and
`F_2` for `Q(sqrt -2)`, a point of the orbit then satisfies

    2P = Q^+ + M,   Q^+ = P + sigma^2 P in E(F_2),   M = P - sigma^2 P in the minus part for F/F_2,

where `Q^+` is one of the finitely many halves of the pinned sum, and the minus part is the
Mordell–Weil group of the twist of E by the element defining `F` over `F_2`: of rank one, proved and
saturated, with generator of height 7.1168638516. So `M` runs over multiples of that generator, and
the exact test decides. The search over `|k| <= 4` leaves one curve, at `k = ±1`, in two
presentations that are the same curve over Q — the isomorphism is exhibited as at X_0^6(187).

**The model is hyperelliptic with no `y^2 = f(x)` over Q**: the conic is `X1^2 + 2X2^2 + 22X3^2 = 0`,
again the one ramified at 2 and at the real place, with

    y^2 = 46521 s^4 + 32184 s^3 z + 101968 s^2 z^2 + 91392 s z^3 + 41984 z^4.

The quotient by the residual involution is `y^2 = A(s, 1)`, whose Jacobian is 1938h1 as predicted;
the counts agree with the trace formula at the 41 good primes below 200 where this model is smooth
and over `F_{p^2}` at 7 primes (the model is singular at 11, which is not a prime of `DN`, an artefact
of the coordinates). **Uniqueness rests on the search window**, since the order-four characters are
not computed; the value found is the smallest nonzero one. This agrees with the Nelson–Shimizu lift,
which reached the same conclusion numerically, including that the canonical conic has no rational
point, by a route sharing no step with this one.

**Checks that were predictions, not inputs.** For every model: the genus is 3, as the genus
formula says; the primes of bad reduction are exactly the primes of `DN` (for the two models over
a conic, the primes of the discriminants of the two equations, apart from the 59 above); the point
counts agree with the trace formula at every good prime below 200 and over `F_{p^2}` at every good
prime below 40, which in genus 3 determines the zeta function; the quotient by the residual
Atkin–Lehner involution is the predicted elliptic curve; the automorphism group has order 4, or 8
at X_0^1326(1) and X_0^210(13) where `W` has index 4 in the Atkin–Lehner group; and the fields of
definition of the fixed points are the fields of the branch places, rational exactly at the places
of degree 1. At X_0^1326(1) the L-polynomial at every good prime below 60 equals the product of
those of 1326b2, 1326c1 and 1326d1, and at X_0^210(13) that of 2730u1, 2730s1 and 210a2. The
models of X_0^39(11) and X_0^55(7) are Q-isomorphic to the curves fitted numerically by the
Nelson–Shimizu lift (campaign branch, `vvdata/weyl-campaign/nelson/hyp16/`), which shares no step
with this construction; at X_0^6(185), where that fit could not recognise two of its invariants,
the exact model agrees with the numerical values to twelve digits.

**⚠ Two of the eighteen are candidates rather than proved.** At X_0^6(385) and X_0^714(5) the
search over the isogeny class is incomplete: Magma does not prove the Mordell–Weil group of the
twist on every curve of the class, even at `Effort := 2`, so a curve with isogenous Jacobian is
not excluded, and the trace formula cannot separate those. Everything else about the two models
checks out. The way to settle them without a new Mordell–Weil computation is to use the isogeny:
if `phi : E -> E'` has degree `m`, the generator of the twist of `E'` is `(c/m)·phi(G)` for some
integer `c` between 1 and `m`, so its height is `c^2/m` times the known one, and the height
comparison on `E'` can be run over those finitely many values.

**What rests on a theorem.** Theorem 1.5 of Cai–Shu–Tian in its second special case, with the
admissibility rule for the characters at the primes dividing both the conductor and the
discriminant, which is that theorem's own test-vector condition: at a ramified prime dividing the
level the involution induced by a uniformiser of `K_v` is the Atkin–Lehner involution `w_v`, and
the theorem asks that the form be a `chi(pi_v)`-eigenvector for it. The Manin constant of the
optimal curve is 1 for squarefree conductor (Česnavičius, *The Manin constant in the semistable
case*, Compositio Math. 154 (2018)), and the LMFDB records it as 1 for every curve used here. The
oldform multiplier above is derived, except at a branch place where the prime is not ramified.
Given those, and given the window at X_0^15(38) and the isogeny-class caveat at the two
candidates, the models are proved.
