# X_0^39(2): a malformed Borcherds form, and the closed form stepping over it

2026-08-29. Found during the Tier 0 backlog sweep, which ran `gtsweep` on the six
bases that used to die inside `RationalNumber`.

## Symptom

`39_2` no longer crashes in `RationalNumber`. It now trips a better guard:

    Runtime error: M0MultiplierExact: class-constancy violated

and, once that guard is made non-fatal, a second one:

    Runtime error: M0MultiplierExact: constant term not real

None of the six Tier 0 bases crashed in `RationalNumber`, so that prediction
holds 4/4 among those that got far enough to say.

## Root cause: form -2's principal part is not integral

`ppdump.m` extracts the infinity-side principal parts without running
`ValuesAtCMPoints`, which is where `39_2` dies. Form -1 is healthy — small
integers throughout. Form -2 is not (`pp_39_2.log`):

    PP -2 78 2
    PP -2 12 -5241430478785720235151540284532191436262054646951/32768
    PP -2 11  8927967484682725495183740059243764757141166954191/65536
    ...

Exact rationals of size ~1e44 with denominators 2^15 and 2^16. Schofer's
formula, and Guo–Yang Thm B with it, assume **c_eta(m) in Z for m <= 0**. Form
-2 violates that hypothesis outright, so it is not a legitimate input,
whatever else is true. Its pole order is also 78, far deeper than anything
else in the panel.

The structure suggests a spurious common term that ought to cancel: m = 12, 10,
4 all carry ±A/32768; m = 11, 8, 2 all carry B/65536; and m = 9, 1 are A/32768
displaced by exactly 65536/32768 = 2. Compare the known "huge odd denominators"
pathology recorded against `21_2`'s split-disc breakage — likely the same
family of failure.

## The guard was right; it is NOT a tolerance artifact

The four guards in `VectorValuedForm.m` use an ABSOLUTE tolerance of 1e-25, and
the obvious hypothesis was that they misfire against large quantities. They do
not. Instrumenting the constancy check to report magnitudes (`diag_39_2.log`)
gives, uniformly across all six violations:

    dev = 1.07e-22    mag = 0.0256    rel = 4.18e-21

`mag` is SMALL. So this is not threshold-vs-magnitude. But `rel = 4e-21` is
~50 digits worse than 80-digit working precision should give: the signature of
**catastrophic cancellation**, consistent with summing 1e44-sized terms down to
0.026. **A relative tolerance would not have fixed this.** The guard caught a
real defect rather than rounding.

(Recorded because I predicted the opposite before measuring, and was wrong.)

## The closed form gets both multipliers anyway

Proposition 9.15 covers this base — omega(D) = 2 for D = 39 = 3·13, and N = 2
is prime — and `mult(f) = (1/2) sum_m c(-m) W(m)` needs only the principal
part. Every pathological coefficient sits at m in {1,2,4,8,9,10,11,12}, and
`W(39,2,m) = 0` at **every one of them**. The support rule annihilates the
garbage:

    form -1:  m=6 (c=-2, W=4/3) and m=13 (c=2, W=4/3)   ->  -8/3 + 8/3 = 0
    form -2:  m=6 (c=-2, W=4/3) and m=78 (c=2, W=4/3)   ->  -8/3 + 8/3 = 0

    mult(-1) = 0,   mult(-2) = 0

Exact rational arithmetic, no complex numbers, no class-constancy, no
tolerance, milliseconds — on a base where the numerical assembly cannot run at
all. The closed form is not merely cheaper here; it is **structurally immune to
the pathology**, because the indices carrying it are exactly the indices where
the local factor vanishes.

## What is NOT established

* No independent ground truth at `39_2` confirms <0,0>: the residual
  measurement needs the pipeline that dies. This is a prediction, not a
  verified value.
* The infinity-only form of the pairing is verified on nine bases / 39 forms
  (`multcheck.py`) but is NOT known to be valid in general — Cor 9.14 pairs
  over all cosets, and §8 measures residual 2.08 for an infinity-only
  functional on the 158-monomial probe space. Promoting the closed form from
  cross-check to production path is gated on closing that.
* Why form -2 is malformed is not diagnosed here. Fixing it, or rejecting it
  and rebuilding, is a separate job.

## Suggested follow-ups

1. Reject non-integral principal parts at construction with a clear error,
   rather than letting them reach `M0MultiplierExact` and surface as an
   inscrutable constancy violation.
2. Cross-check `M0MultiplierExact` against the closed form instead of (or as
   well as) the absolute-tolerance guards: it is scale-free, exact, and turns
   two independent routes into a mutual check.
3. Close the infinity-only question, which now gates a real production use.
