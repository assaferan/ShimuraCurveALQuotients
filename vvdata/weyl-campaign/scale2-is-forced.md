---
name: scale-2-is-forced
description: "The scale 2 is NOT a free parameter: w1+w2 = 1/2 in all nine N>1 fits (forced by w3 = -1/2 and weight-sum 0), and = 1 in all six N=1 fits. So the scale is 1/(1/2) = 2 universally -- no q, s or omega dependence possible. This RETIRES the 'scale 2 vs 2^(omega-1)' question. What remains is the sharper one: why is a_0(E_eis) = 1 at N=1 and 0 at N>1 -- with a concrete suspect, the tracked-coset difference."
metadata: 
  node_type: memory
  type: project
  originSessionId: 6c49569e-411f-4d39-8b37-83e757202459
  modified: 2026-08-26T18:51:59.980Z
---

# The scale 2 is FORCED, not fitted (2026-08-27)

Follows the N=1 deconfounder ([[gross-genus-dictionary]],
[[handoff-2026-08-27]]).  Checked across all 15 banked canonical fits:

| regime | w1 + w2 | total (= a_0(E_eis)) | bases |
|---|---|---|---|
| N > 1 | **1/2**, always | 0 | 9 |
| N = 1 | **1**, always | 1 | 6 |

At N>1 the third weight is w3 = -1/2 always and the total is 0, so the
surviving pair is FORCED to sum to 1/2 -- this was already in the s-law
("w3 = -1/2 ALWAYS; w1 + w2 = 1/2 ALWAYS") but its consequence was not
drawn.  The N=1 fits are the same pair normalized to sum 1.  Hence

    scale = 1 / (w1+w2)|_{N>1} = 1/(1/2) = 2,  universally.

**This retires the "scale 2 vs 2^(omega-1)" question**: the scale cannot
depend on q, s, or omega, because it is fixed by two constants of the
s-law that hold at every base.  It also explains why 210_1 and 330_1
(omega(q)=3) agreed with 10_1/14_1/22_1/26_1 (q prime) -- there was
never any room for them to differ.

## The remaining question, sharpened

Why is a_0(E_eis) = 1 at N=1 but 0 at N>1?

**Concrete suspect: the tracked coset.**  The pipeline tracks rho at a
single isotropic coset, and #iso = 2N-1, so the two regimes track
DIFFERENT things:
* N > 1: `est` = the first NONZERO isotropic coset;
* N = 1: only the zero coset exists, so `est` = i0 (the ZERO coset).

(See the `est` selection in eis32fast.m / eis32k.m -- the N=1 branch was
added for 210_1.)  So the scale 2 and the a_0 jump may be an artifact of
WHICH coset is tracked, not arithmetic about the third Gross slot at all.

**Cheap decisive test (not yet run):** take an N>1 base (e.g. 22_3,
M=132, or any small one), force `est` = i0 and refit.  If the weights
come out as 2x the pair with the third term dropped, the scale 2 is
entirely a tracked-coset normalization and the "missing indefinite
slot" story is a coincidence of the two regimes coinciding.  If the fit
is unchanged or fails, the N=1 story stands and the derivation is
genuinely about the indefinite (D,1) slot.

Do this test BEFORE any Siegel-Weil derivation work -- it decides which
object needs deriving.

## Note for the derivation

Attack via the mass law, which is GAUGE-INVARIANT (|w2/w1| = mass ratio
holds on BOTH supports at the D=2q bases), rather than via the weight
formulas, which depend on the "q = largest ramified prime" convention.
