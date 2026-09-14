# The deficit predictor WORKS — 7/7, and ~100x cheaper than running the base

*2026-09-14. `vvdata/weyl-campaign/deficit.m`, validated on four bases it was never calibrated on.*

## Why this matters

Obstruction is currently discovered by running a base for HOURS and reading the failure. Five
bases were found obstructed today purely by accident, in a routine sweep. And the obstruction is
**not predictable from the shape of `(D,N)`** — every `D = 2p` family contains both builders and
obstructed bases, and `p = 71` is obstructed while `73` builds (see the obstructed-rerun README).

`deficit.m` computes `deficit = Ncols(mat) - Rank(ech_basis * mat)`, which depends only on the
weakly holomorphic basis and the divisor matrix — **not on any divisor choice** — so it skips the
CM points, the field-of-definition work and the 96-triple search entirely.

## The validation

Two bases known OBSTRUCTED and two known BUILT, all even `D`, none used to calibrate the script
(its recorded ground truth was only `38_5 -> 1`, `38_7 -> 0`, `34_3 -> 0`):

    base    truth        deficit at P = 102 / 134 / 190 / 266     total time
    142_1   obstructed        1     1     1     1                    76 s
    158_1   obstructed        1     1     1     1                   146 s
    146_1   BUILT             0     0     0     0                    66 s
    194_1   BUILT             0     0     0     0                   271 s

**4/4 correct; 7/7 including the recorded ground truths.** The deficit is also INVARIANT in the
pole order, which is the script's own criterion for a genuine deficit as against an artifact.

⇒ **~100x cheaper than a pipeline run** (minutes vs hours), and it answers the only question that
matters before committing to a base.

## ⚠⚠ THE LOWEST RUNG LIES, AND IT LIES TOWARD FALSE POSITIVES

At the first pole order (`P ~ 55-73`) the four gave `4, 3, 1, 2` — and **`146_1`, which BUILDS,
reads deficit 1 there**. Reading the first rung alone would condemn a perfectly good base.

    P = 55   146_1  rows 39  cols 19  nds 18  rank 18  deficit 1   <-- FALSE POSITIVE
    P = 102  146_1  rows 86  cols 33  nds 32  rank 33  deficit 0   <-- correct

⇒ **Require `P >= 102` AND require the value to be stable across at least two rungs.** The whole
diagnostic is the invariance, not any single number.

## ⚠ ODD `D`: THE SHORTCUT IS DEAD, AND THE 0-SIDE BLOCK IS GENUINELY REQUIRED (2026-09-14)

A tempting shortcut: the block `deficit.m` omits for odd `D` only ever `VerticalJoin`s **rows** onto
`coeffs_trunc`, and more rows can only raise the rank, so `deficit = Ncols - Rank` can only FALL.
The reported odd-`D` value is therefore an **upper bound** — so a reported **0** would still be
valid (a true deficit cannot be negative), and "clear" is exactly the verdict the screen needs.
No code change required.

**MEASURED, AND IT DOES NOT WORK.** Six odd-`D` bases that all BUILD (true deficit 0):

    15_1 -> 20    39_1 -> 19    51_1 -> 20
    55_1 -> 22    57_1 -> 19    21_2 -> 10

The overestimate is ~20, so **no odd-`D` base will ever read 0** while the block is missing. The
shortcut is useless in practice. (Consistent with the recorded `65_2`: 5, 6, 6, 9.)

⇒ **The 0-side block must actually be implemented.** It is `BorcherdsForms.m:876-978`, 102 lines,
self-contained by its own comment (depends only on `m_idx` plus `D0, n0, nE0, t,
eta_quotients_oo, Xstar`; produces `mat_0_oo` and `relevant_ds_0_oo`).

⚠ **EXTRACT IT, DO NOT COPY IT.** Duplicating 102 lines of submatrix slicing, a kernel solve, two
`coeffs_to_divisor_matrix` calls and an in-place recombination of `ech_etas_0` invites exactly the
drift this repo keeps paying for — and the block's own comment warns those lines "must move
TOGETHER" or you recombine an already-recombined list. Make it a file-local function and `import`
it, as `deficit.m` already does for `basis_of_weakly_holomorphic_forms`.

⚠ **COST:** it refactors the hottest path, inside a memoised loop (`max_pole_order_0`) whose
hoisting was itself a measured optimisation verified by checksumming across 336 triples at `65_2`.
Getting it wrong is a correctness regression on every base, so it needs a FULL SUITE run (~4 h)
before it can be trusted. Deferred deliberately on 2026-09-14; the even-`D` screen (105 targets)
works today and is the larger half.

## Limits

* ⚠ **EVEN `D` ONLY.** For odd `D`, `BorcherdsForms` joins a 0-side block (`coeffs_0_oo`) that this
  script does not compute, so the deficit is an OVERESTIMATE and drifts with the pole order instead
  of staying constant (`65_2`: 5, 6, 6, 9). Do not read odd-`D` numbers until that block is added.
  All four bases here are even `D`.
* Four points plus three recorded ones is not proof. It is, however, enough to justify screening
  the backlog with it rather than attempting bases blind.

## ⇒ The obvious use

**Screen the 180 fresh reachable targets before launching any of them.** At ~1-5 minutes each that
is a few hours of one machine, against days of pipeline runs that end in "Failed to find all
Borcherds forms". It also converts the obstructed class from "whatever we have tripped over" into
something measurable — the current figure, 54, is a lower bound reached entirely by accident.
