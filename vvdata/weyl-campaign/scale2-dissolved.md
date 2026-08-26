# THE SCALE 2 DISSOLVES: it is the TRACKED COSET, not arithmetic

Measured 2026-08-27, immediately after [[scale-2-is-forced]] flagged the
tracked coset as the suspect.  Test as preregistered in that note: force
`est` = i0 (the zero isotropic coset) at an N>1 base and refit.

## The result -- both bases, unambiguous

| base | tracked coset | fit | sum |
|---|---|---|---|
| 15_2 | nonzero iso (default) | -1/2 Gr(5,2) + 1 Gr(5,6) - 1/2 Gr(30,1) | 0 |
| 15_2 | **ZERO** (EST:=0) | **-1 Gr(5,2) + 2 Gr(5,6)** | **1** |
| 22_3 | nonzero iso (default) | -1 Gr(11,3) + 3/2 Gr(11,6) - 1/2 Gr(66,1) | 0 |
| 22_3 | **ZERO** (EST:=0) | **-2 Gr(11,3) + 3 Gr(11,6)** | **1** |

Forcing the zero coset at an N>1 base reproduces the ENTIRE N=1
phenomenology: the third Gross term drops out, the surviving pair
doubles, and the weight sum goes 0 -> 1.  At 22_3 the zero-coset fit is
literally (-2,+3) -- the same numbers as all six N=1 bases.

**So the scale 2 is a normalization convention, not arithmetic.**  There
is nothing to derive.

## Two of this session's readings are thereby REFUTED

1. "The doubling tracks N=1" (the deconfounder's verdict) is not the
   deep story either.  N=1 does not CAUSE anything -- it merely removes
   the choice, because #iso = 2N-1 = 1 leaves the zero coset as the only
   option.  The real variable is which coset is tracked.
2. "The third term is absent at N=1 because (D,1) is indefinite (formal
   mass 0)" is FALSE.  At 15_2 the third slot (30,1) is perfectly
   definite (1 class) and it still vanishes under EST:=0.  The term
   disappears because of the tracked coset, not because the slot is
   indefinite.

## What this does and does not damage

UNAFFECTED (all scale-invariant or normalization-independent):
* The mass-ratio law |w2/w1| = mass2/mass1 -- ratios are scale-invariant.
  This now EXPLAINS why it held 15-for-15 across both regimes and on both
  supports: it is the gauge-invariant content of the fits.
* The SUPPORT verdicts, including the 210_1/330_1 re-split (q = D/2) and
  the exclusion of the competing supports -- scaling all weights does not
  change which columns are used.
* Every residual, rank certificate and gauge check.

NEEDS RESTATING:
* Absolute weights are NOT comparable across the two regimes.  The nine
  N>1 bases are in the nonzero-isotropic-coset normalization; the six
  N=1 bases are in the zero-coset normalization.  The s-law's weight
  formulas apply to the former; the N=1 weights are those same values
  renormalized (doubled, third term dropped).
* The preprint's s-law section must state which coset the weights refer
  to, or quote the mass-ratio form, which is convention-free.

## Provenance / a bug this caught

The control run is what exposed all of this.  A first attempt at the
test had eis32k setting the c=0 (T^k) word's tracked component to 1
unconditionally; that is right only when est = i0, and injects a
spurious unit entry at N>1 (control resid 0.118, |rho| 1.342 vs the
banked 0.894).  Fixed to `(est eq i0) select 1 else 0`; the controls
then reproduce the banked runs exactly (|rho| 0.8944271910 at 15_2,
0.7385489459 at 22_3, both matching to 10 digits) and the 210_1 gate is
unchanged at 1.0087e-42.  **No banked N=1 result was affected** -- there
est = i0 and the old value was correct.
