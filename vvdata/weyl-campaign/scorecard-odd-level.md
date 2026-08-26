# SCORECARD: the odd-level N=1 bases -- measured, same day

Preregistration: `prereg-odd-level.md` (banked 2acaf8a, before any
Gross-span run). Four bases, four verdicts.

## The measurements

All four pass the word-by-word FFT control and the residue pairing:

| base | q, s | \|A\| | CTL worst dev | SOLVE resid | \|rho\| | Eis dim | rank |
|---|---|---|---|---|---|---|---|
| 15_1 | 5, 3 | 450 | 1.32e-79 | 8.00e-42 | 2.1909 | 7 | 7 |
| 21_1 | 7, 3 | 882 | 1.37e-79 | 6.12e-41 | 2.1381 | 7 | 7 |
| 33_1 | 11, 3 | 2178 | 2.51e-79 | 3.24e-42 | 2.0889 | 7 | 7 |
| 35_1 | 7, 5 | 2450 | -- | 9.70e-43 | 2.0284 | 7 | 7 |

Regression gate on the untouched `r2 = 3` path: 15_2, worst
`|closed - FFT| = 7.25e-80` over all 144 words.

## The fits (each restricted support is aliasing dim 0 -- uniquely pinned)

| base | canonical support | weights | predicted | competitor | weights |
|---|---|---|---|---|---|
| 15_1 | (5,1),(5,3) | **(-1, +2)** | (-1, +2) | (3,1),(3,5) | (-1/2, +3/2) |
| 21_1 | (7,1),(7,3) | **(-1, +2)** | (-1, +2) | (3,1),(3,7) | (-1/3, +4/3) |
| 33_1 | (11,1),(11,3) | **(-1, +2)** | (-1, +2) | (3,1),(3,11) | (-1/5, +6/5) |
| 35_1 | (7,1),(7,5) | **(-1/2, +3/2)** | (-1/2, +3/2) | (5,1),(5,7) | (-1/3, +4/3) |

Weight sum is 1 on every fit, canonical and competitor, as the `N=1`
(zero-coset) normalisation requires.

## SCORECARD

**P1 (2x the s-law): HIT, 4 for 4** -- and more than predicted. Every
*competitor* support also obeys `w = (-2/(s'-1), (s'+1)/(s'-1))` read
with ITS own `(q', s')`:

    (3,5) : 2 x (-1/4,  3/4)  = (-1/2, +3/2)   [15_1]
    (3,7) : 2 x (-1/6,  2/3)  = (-1/3, +4/3)   [21_1]
    (3,11): 2 x (-1/10, 3/5)  = (-1/5, +6/5)   [33_1]
    (5,7) : 2 x (-1/6,  2/3)  = (-1/3, +4/3)   [35_1]

exactly. The law is not attached to a choice of `q`; it holds for
either factorisation of a two-prime `D`.

**P2 (constant `(-2,+3)`): REFUTED, decisively.** Not one of the eight
fits above is `(-2,+3)`. The six banked `N=1` bases looked like a
constant only because all six had `s = 2`. The missing cell settles it.

**P3 (the mass-ratio law): HIT, 8 new confirmations.** `|w2/w1|` equals
the mass ratio on every support, canonical and competitor:

    15_1: 2 = (1/6)/(1/12)     |  3 = (1/8)/(1/24)
    21_1: 2 = (1/4)/(1/8)      |  4 = (1/6)/(1/24)
    33_1: 2 = (5/12)/(5/24)    |  6 = (1/4)/(1/24)
    35_1: 3 = (3/8)/(1/8)      |  4 = (1/3)/(1/12)

All masses were computed in advance (`massodd.m`, in the
preregistration). The law now stands at **19 bases**, across both
normalisations, both parities of `DN`, and every support that carries
`E_eis`. It remains the one statement here invariant under all the
conventions.

**P4 (support NOT degenerate at odd D): REFUTED** -- flagged in advance
as the least confident, and rightly. Both ramified primes carry
`E_eis`, each at aliasing dim 0, on all four bases. So degeneracy is
not special to `D = 2q`; it is what two-prime discriminants do. This
sharpens rather than weakens the four-prime result: at `210_1` the
competing support was genuinely EXCLUDED (rank 87 vs 88), and at
`330_1` likewise (135 vs 136). **Two primes cannot select; four primes
do.** That is now a clean statement about where the re-split rule has
content -- and it is the reason `1155_1` is the test that matters.

**P5 (level constants recur at 15_1): REFUTED as stated, my error.**
I conflated the eta-monomial level `M = 2DN` with the form level
`4DN`: 15_1 and 15_2 share `M = 60` but their weight-3/2 spaces are at
levels 60 and 120. The odd-level family has its own uniform pair,
**Eisenstein dimension 7 and constants-rank 7 at all four bases**
(levels 60, 84, 132, 140), with aliasing 1 on the full pool and 0 on
either restricted support. So the "base-independent constants" are per
level-family, as the corrected preprint already says.

## What this settles

The `N=1` doubling is uniformly `2 x (s-law)`. Together with the
`EST:=0` experiment (forcing the zero coset at an `N>1` base
reproduces the `N=1` shape), the two regimes are now one law in two
normalisations, with nothing left to derive on the scale.
