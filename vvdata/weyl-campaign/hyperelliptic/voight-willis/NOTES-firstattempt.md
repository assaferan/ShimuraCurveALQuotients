# Voight–Willis feasibility spike — 2026-10-15

Go/no-go on computing holomorphic differentials of Shimura curves analytically, to feed Max
Noether's count (FilterByCanonicalQuadrics) and decide the 34 open quotients.

**The paper itself was not consulted** (no network access to it from here). The method was
reconstructed from its description, so nothing below is a check against the source.

## What Magma already provides — all of it works

`FuchsianGroup`, `FundamentalDomain`, `Group` (presentation), `WordProblem`, `ShimuraReduceUnit`,
`CMPoints`, `Signature`, `FuchsianMatrixRepresentation`, `UnitDisc`, `PlaneToDisc`/`DiscToPlane`,
and the undocumented attributes `G`ShimFDGens` (side-pairing generators as quaternions) and
`G`ShimFDDisc` (domain vertices in the disc). Group elements act on disc points with `*`.
`SingularValueDecomposition` and `NumericalRank` exist (SVD only over the reals, so a complex
matrix must be expanded to the real 2m x 2n form — done in `vwlib.m`).

⚠ `jParameter` requires `G`TriangleBool` and is therefore useless here: **X_0^6(1) has signature
(0; 2,2,3,3), four elliptic points, so it is NOT a triangle group.** The directive's rung-1 plan
of recovering a Hauptmodul through `jParameter` cannot be run.

## The expansion layer I built (`vwlib.m`)

Disc matrix = Cayley conjugate of the upper-half-plane matrix, normalised to determinant 1; this
comes out **exactly in SU(1,1)** (checked: the ratios W[2][2]/conj(W[1][1]) and
W[2][1]/conj(W[1][2]) are 1 to full precision), so the automorphy factor j(g,w) = cw+d is
unambiguous for even weight. Impose f(gw) = j(g,w)^k f(w), truncating f = sum_{n<N} a_n w^n, and
read the space off the numerical nullspace.

**Sampling is the crux.** Sample on each generator's isometric circle |cw+d| = 1 *restricted to
the arc inside the domain radius R* (keep w with |w| <= R and |gw| <= R). Sampling whole isometric
circles instead reaches |w| = 0.966 at D = 6, where 40 terms give a truncation error of 0.25 and
the method returns pure noise. With the restriction, max|w| drops to 0.552 and the error is R^N.

**Column normalisation is required.** Untreated, the trailing columns have entries of order R^n
and are numerically negligible, producing spurious near-null directions whose *count grows with
the truncation length* (2 at N = 40, 4 at N = 80). `NormaliseColumns` removes that artefact.

## Results

Dimensions of M_k for X_0^6(1), signature (0; 2,2,3,3), against
dim = (k-1)(g-1) + sum_i floor((k/2)(1-1/e_i)), at precision 80 and N = 60 terms (R = 0.554):

| k | expected | got | smallest singular values |
|---|---|---|---|
| 4 | 1 | 4 | 1.3e-82, 1.4e-78, 1.8e-62, 6.0e-62, then 7.3e-6 |
| 6 | 1 | **1 OK** | 4.8e-16 (= R^N), then 2.0e-5 |
| 8 | 1 | 4 | four at machine zero, then 1.9e-6 |
| 10 | 1 | **1 OK** | 2.8e-15, then 1.8e-5 |
| 12 | 3 | 4 | four at machine zero, then 5.2e-6 |
| 14 | 1 | **1 OK** | 4.9e-15, then 1.2e-5 |

So **k = 2 mod 4 is correct**, with the vanishing singular value sitting exactly at the predicted
truncation level R^N — a strong self-consistency check — while **k = 0 mod 4 carries four spurious
*exact* null vectors** (machine zero, not truncation level). Unexplained. Not column scaling, not
truncation length, and not the automorphy convention: six variants of j (j, -j, i*j, 1/j, conj j,
j/|j|) all give the same pattern, and all six side-pairing generators do contribute sample points
(7, 38, 38, 17, 38, 38 of 400).

**Weight 2 — the weight we actually need, where dim = genus.** Precision 80, N = 60:

| curve | genus | R | R^N | got | verdict |
|---|---|---|---|---|---|
| X_0^6(1) | 0 | 0.554 | 4.2e-16 | 1 | MISMATCH |
| X_0^14(1) | 1 | 0.707 | 8.8e-10 | 1 | number right, but **no gap** (ratio 1.38; smallest sv 1.6e-4, far above truncation) — not a pass |
| X_0^15(1) | 1 | 0.805 | 2.2e-6 | 1 | same, gap 3.7, smallest sv 2.4e-5 |
| X_0^26(1) | 2 | 0.867 | 2.0e-4 | 1 | MISMATCH |

## Why it does not scale — two independent blockers

1. **The domain radius R rises toward 1 with the curve**: 0.554, 0.707, 0.805, 0.867 at
   D = 6, 14, 15, 26. Terms needed go like log(eps)/log(R), so at D = 26 sixty terms already leave
   only 2e-4 of resolution and the answer is wrong. Our targets are two orders of magnitude larger
   in level. The cost also shows in the SVD: 7 s, 42 s, 165 s, 260 s for those four curves.
2. **The prerequisite itself stops first.** `FundamentalDomain` takes 248 s at D*N = 210 (102
   sides) and **did not finish in 30 minutes at D*N = 770**, which is the *smallest* of the 34 open
   curves by level (X_0^77(10)); the largest is D*N = 4130. This caps the approach before any
   expansion layer matters.

On top of both, the Atkin–Lehner layer (extending the group by the involutions to get X/W rather
than X) is extra work not in the published method, and every verdict would be numerical.

## Verdict: NO-GO

Scripts: `vwlib.m` (layer), `explore*.m`, `dbg*.m` (diagnostics), `rung1.m`, `rung1b.m`, `wt2.m`,
`time.m`/`time.log` (prerequisite timings).
