# Power series expansions of modular forms on Shimura curves (Voight–Willis)

A second, independent way to produce **equations for Shimura curves and their Atkin–Lehner
quotients**, alongside the Guo–Yang / Borcherds route the rest of this repository uses.

Reference: J. Voight and J. Willis, *Computing power series expansions of modular forms*, in
Computations with Modular Forms (Springer, 2014); author's copy at
`jvoight.github.io/articles/ModFormPaper-090317.pdf`. The implementation follows its Section 3.

Where the Borcherds route builds a model from an eta-quotient basis obtained by enumerating lattice
points of a polytope, this one expands modular forms numerically in a disc coordinate on the
Fuchsian group's fundamental domain, reads off a basis of forms of each weight, and builds the model
from their ratios. The two are limited by completely different things, which is the point of having
both:

| | Borcherds / Guo–Yang | Voight–Willis |
|---|---|---|
| what blocks it | the number of divisors of M = 2DN (the Normaliz wall at 24), and the m = 0 multiplier when N > 1 | the size of the fundamental domain, which grows with D·N |
| exactness | exact equations | numerical, so a verdict needs an error statement or an exact post-check |
| conditionality at N > 1 | rests on the m = 0 term | none |

So a base can be out of reach of one and inside the other. The bases behind the divisor wall sit at
D·N = 330, 462, 510 and upward with many divisors, while this method does not care how many
divisors D·N has. That is the opportunity to test.

## Status

**Validated at X_0^6(1)** (signature (0; 2,2,3,3), domain radius 0.5435):

* `dim M_k` for k = 2, 4, …, 14 reproduces the orbifold dimension formula (0, 1, 1, 1, 1, 3, 1),
  with the vanishing singular values separated from the rest by 1e6 to 1e9;
* the Hauptmodul reproduces Elkies (*Shimura curve computations*, arXiv:math/0005160, Sec. 3.1 and
  Table 1, transcribed with its corrections in `tests/CMFieldsElkies.m` on `main`): the cross-ratio
  of the elliptic points is −1 to 7.7e-15, and his t = 0 at discriminant −24 comes out as 1e-34;
* a check **not** forced by symmetry: the valence formula gives deg div(f_4) = 2/3, carriable only
  by the two order-3 points at 1/3 each, and the computed weight-4 form vanishes at exactly those
  two, at relative magnitude 8.4e-21 against 0.503 and 1.000 elsewhere.

## What has to happen next, in order

1. **The Atkin–Lehner layer.** Everything so far is on X_0^D(N) itself. Every model we store, and
   every curve we care about, is a quotient X/W, so the group must be extended by the Atkin–Lehner
   elements and the forms cut down to the W-invariant ones. This is beyond the published method and
   is the first real piece of work.
2. **Agreement with the models we already have.** The decisive test is not a dimension count: build
   models from the expansions at bases where `data/models/models_D_N.m` already holds equations, and
   check they agree with the stored ones. Note the model files are scattered across branches, so
   assemble a merged directory first (see `../hyperelliptic/open-cases/have.m`). Compare through the
   checks the repository already trusts: the genus from the genus formula, point counts over F_p and
   F_{p^2} against the trace formula (`ComputePointsViaTrace`, as `tests/CRVTraceCounts.m` does),
   and curve isomorphism by construction. Agreement across many bases is what would make this a
   second source rather than a second guess.
3. **Beyond the Borcherds reach.** Only then try the bases the polytope cannot do — the 24-divisor
   levels — and the quotients still undecided for hyperellipticity.

## The wall, measured

`FundamentalDomain` is the prerequisite and the bottleneck:

| D·N | result |
|---|---|
| 210 | 248 s, 102 sides, presentation 18 generators |
| 770 | not finished in 30 min, on both the Mac and lovelace |

A 12 h run at 770 was launched on 2026-10-15 (`lovelace:~/gymodels/vwdom770.log`) to replace
"timed out" with a number. Note the comparison that matters is against D·N, not against the
polytope level M = 2DN: the first bases behind the divisor wall are at D·N = 330, 462 and 510,
which are between the two rows above, so whether they are reachable is genuinely open.

## Traps, each of which cost a day

* **Magma's `SingularValueDecomposition(M)`** on an m × n matrix returns (m × n, m × m, n × n), and
  the singular values are the diagonal of the **first** return, descending; the second is the left
  orthogonal matrix. Reading its diagonal instead gives numbers bounded by 1 and of arbitrary sign,
  and produced an entirely spurious earlier conclusion that the method does not work.
* **Equilibrate the columns.** They scale like rho^n, so without it every singular value sits at
  rho^N and no gap is visible at all.
* **Magma's own `CMPoints` is broken at D = 6**, failing with a `RingClassGroup` coercion error
  inside `package/Geometry/GrpPSL2/GrpPSL2Shim/reciprocity.m` line 29. Take instead the fixed point
  of a trace-zero element of reduced norm m.
* The undocumented `ShimFDGens` and `ShimFDDisc` hold the side-pairing generators and the domain
  vertices; the Cayley conjugate of the half-plane matrix lands exactly in SU(1,1), so the
  automorphy factor is unambiguous at even weight.

## Files

`vw3.m` is the working library; `dims.m` computes dimensions of M_k, `haupt.m` builds a Hauptmodul
as a ratio of weight-12 forms, `cmchk.m` evaluates it at a CM point, `predict.m` is the
valence-formula prediction. `vw2.m`, `run2.m`, `shape*.m` and `dbg_*.m` are development history,
kept only because they record what was tried. `vwlib.m` and `NOTES-firstattempt.md` are from the
first attempt and are superseded, except that `time2.m` still serves for timing the prerequisite.
`NOTES.md` is the working log.
