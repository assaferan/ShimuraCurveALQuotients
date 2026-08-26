---
name: scale-2-is-forced
description: "RESOLVED 2026-08-27: the scale 2 DISSOLVED -- it is the tracked isotropic coset, a normalization convention, not arithmetic. Forcing est=i0 at an N>1 base reproduces the whole N=1 phenomenology (third term drops, pair doubles, sum 0->1; 22_3 gives literally (-2,+3)). Refutes BOTH 'N=1 causes it' and 'the third slot is absent because indefinite'. Mass-ratio law and all support verdicts UNAFFECTED (scale-invariant); absolute weights are NOT comparable across regimes."
metadata: 
  node_type: memory
  type: project
  originSessionId: 6c49569e-411f-4d39-8b37-83e757202459
  modified: 2026-08-26T19:17:00.647Z
---

# The scale 2: forced, then DISSOLVED (2026-08-27, banked 9c3153a + 5385deb)

## Step 1 — it was never a free parameter

Across all 15 banked canonical fits: at N>1 the pair sums to **1/2**
always (forced by w3 = −1/2 and weight-sum 0, both s-law constants); at
N=1 it sums to **1**. So scale = 1/(1/2) = 2 identically, with no room
for q/s/ω dependence. That alone retired the "scale 2 vs 2^(ω−1)"
question and explained why composite and prime N=1 bases had to agree.

## Step 2 — THE ANSWER: it is the tracked coset

The pipeline tracks ρ at ONE isotropic coset and #iso = 2N−1, so the
regimes track different objects: N>1 uses the first NONZERO isotropic
coset, N=1 is forced to the ZERO coset. Test (`EST:=0` override in
eis32k.m) — force the zero coset at an N>1 base:

| base | tracked coset | fit | sum |
|---|---|---|---|
| 15_2 | nonzero (default) | −1/2 Gr(5,2) + 1 Gr(5,6) − 1/2 Gr(30,1) | 0 |
| 15_2 | **zero** | **−1 Gr(5,2) + 2 Gr(5,6)** | **1** |
| 22_3 | nonzero (default) | −1 Gr(11,3) + 3/2 Gr(11,6) − 1/2 Gr(66,1) | 0 |
| 22_3 | **zero** | **−2 Gr(11,3) + 3 Gr(11,6)** | **1** |

The zero coset reproduces the ENTIRE N=1 phenomenology at an N>1 base —
at 22_3 literally (−2,+3), the same numbers as all six N=1 bases.
**The scale 2 is a normalization convention. There is nothing to derive.**

## What this refutes (both are mine, from earlier the same day)

* "The doubling tracks N=1" ([[gross-genus-dictionary]]'s deconfounder
  verdict): N=1 does not cause it — it only removes the choice, since
  #iso = 1 leaves the zero coset as the only option.
* "The third term is absent at N=1 because (D,1) is indefinite, formal
  mass 0": FALSE. At 15_2 the third slot (30,1) is definite (1 class)
  and still vanishes under EST:=0.

## What survives untouched

* **The mass-ratio law** |w2/w1| = mass2/mass1 — ratios are
  scale-invariant. This now EXPLAINS its 15-for-15 run across both
  regimes and on both supports: it is the gauge-invariant content.
* **Every support verdict**, including the 210_1/330_1 re-split (q = D/2)
  and the rank-exclusions of competing supports — rescaling all weights
  does not change which columns are used.
* All residuals, rank certificates, gauge checks.

## What must be restated

Absolute weights are NOT comparable across regimes: the nine N>1 bases
are in the nonzero-isotropic normalization, the six N=1 bases in the
zero-coset one. The s-law weight formulas apply to the former. **The
preprint must either say which coset the weights refer to, or quote the
convention-free mass-ratio form.**

## Bug the control caught (why controls matter here)

eis32k set the c=0 (T^k) word's tracked component to 1 unconditionally —
right only at est = i0, a spurious unit entry at N>1 (control resid
0.118 vs banked ~5e-42). Fixed to `(est eq i0) select 1 else 0`;
controls then reproduce banked runs exactly (|ρ| 0.8944271910 at 15_2,
0.7385489459 at 22_3) and the 210_1 gate is unchanged. **No banked N=1
result was affected** — there est = i0 and the old value was correct.

Related: [[gross-genus-dictionary]], [[eis32k-closed-form-driver]],
[[handoff-2026-08-27]].
