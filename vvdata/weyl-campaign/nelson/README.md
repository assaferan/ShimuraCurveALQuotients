# Evaluating modular forms on Shimura curves by the Shimizu lift (Nelson)

A third way to get at a modular form on a Shimura curve, alongside the Borcherds / Guo–Yang
eta-quotient route the repository uses for equations and the Voight–Willis power series in
`../voight-willis/`.

Reference: P. Nelson, *Evaluating modular forms on Shimura curves*, Math. Comp. **84** (2015)
2471–2503, arXiv:1210.1243. What is implemented is Theorem 3.1 of that paper and nothing else.

Voight and Willis cite Nelson and write that it would be interesting to compare the two
techniques. This directory is that comparison, and the short answer is that the two agree
wherever both can be run, and that Nelson's cost is governed by `D·N` while Voight–Willis' is
governed by the circumradius of the fundamental domain — which is why one of them reaches the
bases of our open quotients and the other does not.

## What the method computes, and what it needs

For a Hecke eigenform `f_B` of even weight `k` on `X_0^D(N)`, specified by the eigenvalues of its
Jacquet–Langlands correspondent `f`, a newform of level exactly `DN`, Theorem 3.1 says

    (y1 y2)^{k/2} conj(f_B(z1)) f_B(z2)
      = sum_{d | DN} c_f(d) sum_{alpha} a_f(d nrd(alpha)) e^{i k arg X(alpha)} W_k(4 pi |X(alpha)|)

where the inner sum runs over the elements of positive reduced norm of a lattice `R^(d)` between
an Eichler order `R` of level `N` and its dual, each `alpha` is carried into `M_2(R)` by
`sigma_{z1}^{-1} iota(.) sigma_{z2}`, and

    X([a b; c d]) = (a + d + i(b - c))/2,      P([a b; c d]) = (a^2+b^2+c^2+d^2)/2,
    c_f(d)        = mu((d,N)) / (d a_f(d)),
    W_k(x)        = 2^{1-k} sum_{n>=1} n^k ( n x K_{k-1}(nx) - K_k(nx) ).

The `a_f(n)` are the **classical** coefficients, `f = sum a_f(n) e(nz)` with `a_f(1) = 1`.

So the input is: a quaternion algebra of discriminant `D`, an Eichler order of level `N`, a real
splitting, and the Fourier coefficients of one newform. The output is the value of one eigenform
at a time — not a space, which is the structural difference from Voight–Willis.

**An individual value `f_B(z)` is only defined up to a constant of modulus one** (Nelson's
Remark 2.5); it is the product of a value with a conjugate value that the Hecke data determines.
For building a basis of differentials that costs nothing: fixing `z1` once and reading
`conj(f_B(z1)) f_B(z_j)` across `j` gives each form normalised so that its value at `z1` is 1, and
rescaling the members of a basis is a change of coordinate on the canonical image.

**What makes the sum finite** is the identity `2|X(alpha)|^2 = P(alpha) + det(alpha)`. Since `W_k`
decays like `e^{-x}`, truncating at `4 pi |X| <= cut` is the positive definite condition
`P(alpha) <= 2 (cut/4 pi)^2`, so each term of the divisor sum is an enumeration of short vectors
of a rank-four positive definite form. **There is no fundamental domain, no sampling, no
reduction and no matrix anywhere in the method.**

`R^(d)` is built as `R + (DN/d) R*`, which is right because `R*/R` has exponent dividing `DN`;
the code checks `R ⊆ R^(d) ⊆ R*` with index `d^2`, and that the reduced norms on `R^(d)` have
denominators dividing `d`.

## How to run it

All scripts take `name:=value` arguments and must be run with stdin redirected.

    cd worktrees/rerun
    NB=.../worktrees/campaign/vvdata/weyl-campaign/nelson

    # stage 1: agreement with Voight-Willis at a base where both work
    magma -b DISC:=6  K:=4 PR:=30 CUT:=50 TGT:=20 NPTS:=8 $NB/stage1.m   < /dev/null
    magma -b DISC:=15 K:=4 PR:=30 CUT:=50 TGT:=24 NPTS:=9 $NB/stage1.m   < /dev/null

    # negative controls, ONE PER PROCESS (see below for why)
    for m in 0 1 2 4 5; do
      magma -b DISC:=6 K:=4 PR:=30 CUT:=40 TGT:=20 NPTS:=8 MODE:=$m $NB/control2.m < /dev/null
    done

    # stage 2: cost.  cost.m rebuilds the W_k table per base; cost2.m builds it once and
    # enumerates by process, which is what the large bases need.
    magma -b BASES:=6_1,10_1,15_1 K:=2 PR:=30 CUT:=40 CUT0:=32 $NB/cost.m  < /dev/null
    magma -b BASES:=330_1,462_1,510_1 K:=2 PR:=30 CUT:=32 CUT0:=26 NF:=3 $NB/cost2.m < /dev/null

    # stage 3: Noether's count of canonical quadrics, differentials supplied by Nelson
    magma -b DISC:=35 WS:=1 PR:=30 CUT:=36 NPTS:=14 $NB/stage3.m < /dev/null

    # where the eigenvalues have to come from, timed
    magma -b LEVELS:=282,510,1122 CUT:=32 $NB/eig.m < /dev/null

    # candidates for a stage-3 control: genus-3 quotients of X_0^D(1), and among those the ones
    # with no genus-0 quotient above them, which are the non-hyperelliptic ones
    magma -b DMAX:=250 $NB/scan3.m < /dev/null
    magma -b $NB/scan4.m < /dev/null

    # does Theorem 3.1 reach the whole space of differentials at a base with N > 1?
    magma -b "BASES:=39_11:1,11,39,429;141_2:1,3,94,282" $NB/nold.m < /dev/null

    # stage 3 at N > 1 (only the newforms of level exactly D*N; stops if they do not span)
    magma -b DISC:=55 LEV:=7 WS:=1,35,55,77 PR:=30 CUT:=46 NPTS:=18 $NB/stage3n.m < /dev/null

`cmpvw.m` is the earlier ratio form of the stage-1 test and `scan3.m` the search that found the
control pair. `t0.m`–`t3.m` and `twk.m` are the probes the traps below were found with.

## Stage 1: agreement with Voight–Willis

For each eigenform, Theorem 3.1 is evaluated at `m` points and the vector of values is asked to
lie in the span of the Voight–Willis basis of the same space at the same points; with a space of
dimension `g` that is `m - g` independent predictions. The residual is the smallest singular
value of the normalised `(g+1) x m` matrix and the next singular value is its scale.

| base | weight | dim | eigenform | points | cutoff | residual | next sv | gap | predictions |
|---|---|---|---|---|---|---|---|---|---|
| X_0^6(1) | 4 | 1 | the one | 8 | 35 | 2.5e-10 | 1.41 | 5.7e9 | 7 |
| X_0^6(1) | 4 | 1 | the one | 8 | 50 | **1.7e-16** | 1.41 | 8.2e15 | 7 |
| X_0^15(1) | 4 | 2 | first | 9 | 50 | **1.7e-16** | 0.90 | 5.2e15 | 7 |
| X_0^15(1) | 4 | 2 | second | 9 | 50 | **1.4e-16** | 0.73 | 5.0e15 | 7 |

Raising the cutoff from 35 to 50 moves the residual from 2.5e-10 to 1.7e-16, which is what says
the disagreement is the truncation of Nelson's sum and not a systematic difference between the
two methods.

**`X_0^15(1)` is the case with content.** The space is two-dimensional there, so the eigenvalues
are what distinguishes the two forms, and each of the two newforms of `S_4^{new}(Gamma_0(15))`
comes out as its own vector lying in the Voight–Willis plane. At a one-dimensional space the
answer is forced up to a scalar whatever eigenvalues are supplied, so `D = 6` alone would not
have tested the Jacquet–Langlands input at all.

### Checks that do not use Voight–Willis

* **The kernel is rank one.** `|S(z1,z2)|^2 = S(z1,z1) S(z2,z2)` to 1.1e-10 at cutoff 35, that is
  to the truncation level: Cauchy–Schwarz with *equality*. This is the one internal check not
  forced by a symmetry of the truncated sum, and it tests the lattices `R^(d)`, the coefficients
  `c_f(d)`, the values `a_f(n)` and `W_k` all at once, since a wrong ingredient gives a kernel of
  rank greater than one.
* **The diagonal is real and positive**, as `y^k |f_B(z)|^2` must be: imaginary part 1e-44
  against a real part of 7.5e-4, at 40 digits.
* **The transformation law** `S(z1, g z2)/S(z1, z2) = (c z2 + d)^k / |c z2 + d|^k` holds to 1e-29
  for all three generators of the unit group. ⚠ **This one is exact by symmetry of the truncated
  sum** and is not evidence about the truncation: `R^(d) g = R^(d)` for a unit `g`, and
  `sigma_{gz}` differs from `g sigma_z` by a rotation, under which `|X|` is unchanged and
  `arg X` shifts uniformly. It does test the splitting and the two-sidedness of `R^(d)`.

### What the controls control for, and one that does not work

Each is a single input spoiled, run in its own process, at `X_0^6(1)` weight 4 with cutoff 40.
Residual against the Voight–Willis line, correct input first:

| what is spoiled | residual | factor over correct |
|---|---|---|
| nothing | 8.5e-13 | — |
| `a_f(5)` negated (smallest prime not dividing `DN`) | 1.0e-3 | 1.2e9 |
| the newform of level 10 substituted for the one of level 6 | 5.4e-2 | 6.4e10 |
| the divisor sum cut to `d = 1` | 1.9e-1 | 2.2e11 |
| the divisor sum missing only its top term `d = DN` | 2.2e-3 | 2.5e9 |
| every `R^(d)` replaced by `R` | 1.9e-1 | 2.2e11 — **but see below** |

⚠ **The last one is not an independent control.** At a ramified prime `a_f(dt) = a_f(d) a_f(t)`,
so with every lattice equal to `R` the divisor sum collapses to `(sum_d c_f(d) a_f(d))` times the
`d = 1` term, which at `D = 6`, `k = 4` is exactly 2; the computed values are twice those of the
`d = 1` control to every one of their 30 digits. It fires, but it is the same wrong object as the
row above it, scaled. The control that spoils the divisor sum without collapsing it is the
`d = DN` row.

Permuting the lattices (using `R^(DN/d)` where `R^(d)` belongs) cannot be used as a control at
all, because `d nrd(alpha)` then fails to be an integer and the run stops on that assertion —
which is itself a check that the `R^(d)` are the lattices the formula wants.

## Stage 2: what an evaluation costs, and how it grows

The whole expense is the lattice enumeration, and

    number of lattice points  ~  10 * (cut / 4 pi)^4 * sigma_2(DN) / DN,   sigma_2(DN) = sum_{d | DN} d^2

within about 5 percent at every base measured, `DN` from 6 to 4830 and cutoffs 32, 35, 40, 45 and
50. Two further laws, both measured:

* the coefficients `a_f(n)` are needed up to exactly `DN * (cut/4 pi)^2` (predicted and measured
  agree to within 10 at every base);
* **the cutoff buys about 0.4 decimal digits per unit**, with a scatter of about one and a half
  digits between bases and between choices of point, and **no trend in `DN`**: the relative
  change between cutoffs 26 and 32 is 4.5e-8 at `DN = 770` and 1.8e-5 at `DN = 1326`, with
  `DN = 85` (7.0e-8) and `DN = 510` (3.0e-6) in between in neither order.

**Several eigenforms cost one enumeration, not several.** Everything but `a_f(n)` is shared, so a
basis of `g` differentials at one pair of points is one enumeration: three forms take the same
time as one, to within the measurement. That is the figure that decides what a model costs.

Throughput is about 13 800 lattice points per second at 30 digits on one core, flat across every
base measured.

### Measured (weight 2, 30 digits, cutoff 32, three forms per enumeration, one core)

All coefficient rows below are **genuine**: an elliptic curve of conductor `DN` is a
Jacquet–Langlands correspondent and its `a_p` are point counts, so the convergence column means
what it says. Rows using surrogate coefficients are marked and measure **timing and counts only**.

| base | volume | DN | sigma_2/DN | lattice points | seconds | n_max | rel. change 26→32 |
|---|---|---|---|---|---|---|---|
| X_0^85(1) | 10.67 | 85 | 88.7 | 36 741 | 2.4 | 549 | 7.0e-8 |
| X_0^91(1) | 12.0 | 91 | 93.4 | 38 691 | 2.6 | 587 | 2.5e-6 |
| X_0^123(1) | 13.3 | 123 | 136.8 | 56 683 | 3.9 | 796 | 6.9e-8 |
| X_0^142(1) | 11.67 | 142 | 177.5 | 73 234 | 4.9 | 903 | 5.3e-6 |
| X_0^6(35) | 16.0 | 210 | 309.5 | 128 359 | 9.3 | 1 357 | 7.5e-8 |
| **X_0^141(2)** | **46** | 282 | 391.8 | 161 920 | 11.1 | 1 825 | 4.9e-7 |
| X_0^330(1) | 13.3 | 330 | 480.6 | 199 328 | 13.5 | 2 131 | 9.4e-8 |
| **X_0^55(7)** | **53.3** | 385 | 412.0 | 170 891 | 12.5 | 2 494 | 7.0e-8 |
| **X_0^39(11)** | **48** | 429 | 483.5 | 200 507 | 14.1 | 2 780 | 1.8e-7 |
| X_0^462(1) | 20.0 | 462 | 660.2 | 272 787 | 18.8 | 2 994 | 2.3e-7 |
| X_0^510(1) | 21.3 | 510 | 739.2 | 306 494 | 20.0 | 3 302 | 3.0e-6 |
| **X_0^77(10)** | **180** | 770 | 1 029.9 | 427 083 | 32.1 | 4 990 | 4.5e-8 |
| **X_0^6(187)** | **72** | 1 122 | 1 576.7 | 654 094 | 50.4 | 7 262 | 1.7e-7 |
| **X_0^1326(1)** | **64** | 1 326 | 1 859.0 | 771 407 | 53.2 | 8 587 | 1.8e-5 |
| **X_0^6(649)** | **240** | 3 894 | 5 454.6 | 2 262 332 | 182 | 25 247 | 6.5e-8 |
| **X_0^690(7)** | **234.7** | 4 830 | 7 132.5 | 2 958 039 | 232 | 31 298 | 5.0e-8 |

Bases of the 34 open quotients are in bold. At a cutoff of 40 rather than 32 (about three and a
half digits more) every lattice-point count is 2.4 times larger; at 50 it is 6 times larger.

An earlier sweep at cutoff 40 over `DN` = 6, 10, 15, 22, 26, 35, 51, 58 gave the same law with the
same constant and the same precision behaviour, with genuine coefficients at `DN` = 15, 26, 35,
51 and 58 and surrogates at 6, 10 and 22 (`S_2` is zero at those levels, so no genuine weight-2
correspondent exists there and those three rows are timing only).

### Is the cost linear in D·N? Yes, and here is why

For squarefree `DN`, `sigma_2(DN)/DN = prod_{p | DN} (p + 1/p) = DN * prod_{p | DN} (1 + 1/p^2)`,
and that last product lies between 1 and about 1.5. So the law above **is** linear in `DN`, with
a coefficient `10 (cut/4pi)^4 prod (1 + p^{-2})` that is larger when `DN` has small prime factors.
At a cutoff of 40, which buys about ten digits, that is 1030 to 1500 lattice points per unit of
`DN` — the figure of roughly 1200 per unit is right, and the spread around it is the prime
factorisation, not noise.

Voight–Willis, by contrast, is governed by the circumradius `rho` of the fundamental domain
through `NT ~ 1/(-log rho)` and a dense decomposition going as `NT^3`, which is 68 hours and
984 GB already at volume 21. Nelson does not compute a domain at all, and `rho` never enters.
X_0^510(1), the base where that estimate was made, is 306 494 lattice points and **20 seconds**
here.

### What breaks first, and where this stops

I expected the Fourier coefficients to be the wall and **measured that they are not** (`eig.m`,
weight 2, one core): decomposing `S_2(Gamma_0(L))` into newforms and expanding one of them to the
`n_max` the cutoff demands costs

| level | dim S_2 | n_max | newform decomposition | expansion | orbits |
|---|---|---|---|---|---|
| 282 | 45 | 1 839 | 0.20 s | 0.43 s | 5 |
| 385 | 45 | 2 507 | 0.16 s | 0.77 s | 8 |
| 429 | 53 | 2 792 | 0.16 s | 0.97 s | 8 |
| 510 | 101 | 3 318 | 0.25 s | 1.5 s | 8 |
| 770 | 137 | 5 004 | 0.44 s | 3.3 s | 13 |
| 1 122 | 209 | 7 286 | 0.82 s | 7.3 s | 19 |
| 1 326 | 245 | 8 609 | 0.75 s | 10.1 s | 14 |
| 2 002 | 329 | 12 993 | 1.3 s | 23.7 s | 19 |
| 2 730 | 657 | 17 713 | 5.0 s | 48.9 s | 38 |
| 3 894 | 713 | 25 261 | 4.2 s | 95.8 s | 38 |
| 4 830 | 1 137 | 31 331 | 33 s | 186 s | 57 |

So at `DN = 4830`, the largest level in the open range, the eigenvalues cost four minutes and the
evaluation four. Nothing in the method breaks anywhere in the range of the open bases.

Of the three candidate bottlenecks:

* **the enumeration** is linear in `DN` and in hand — 232 s for one evaluation of three forms at
  the largest open base;
* **cancellation** does not happen: the largest single term is 0.3 to 0.5 of the answer, and the
  relative accuracy at a fixed cutoff shows no trend over a factor of 50 in `DN`;
* **the eigenvalues** grow as the table above, the expansion roughly quadratically in `DN`.

Extrapolating those two laws — linear for the enumeration, between quadratic and cubic for the
newform work (the decomposition went from 4.2 s to 33 s between `DN` = 3894 and 4830, a factor of
8 for a factor of 1.24 in level, so that one is accelerating sharply) — the first to become
uncomfortable is the newform decomposition, not far past where we need it. At `DN` of ten to
twenty thousand the decomposition is the problem and the lattice sum is still minutes. That
extrapolation rests on a clean fit for the enumeration (five percent over a factor of 800 in
`DN`) and a poor one for the decomposition — two points, and the second one eight times the
first — so it should be read as "there is headroom above 4830, but not an order of magnitude of
it", not as a number.

**The real limitation is not cost at all; it is the hypothesis of Theorem 3.1.** The theorem is
stated for a newform of level exactly `d_B N`, and at `N > 1` the space of differentials of
`X_0^D(N)/W` also contains oldform copies of newforms of level `DM` for proper `M | N`. How much
of the space the theorem reaches has to be checked base by base, which is what `nold.m` does:

| base | W | genus | W-fixed newforms of level exactly DN | |
|---|---|---|---|---|
| X_0^39(11) | {1,11,39,429} | 3 | 3 | **complete** |
| X_0^55(7) | {1,35,55,77} | 3 | 3 | **complete** |
| X_0^1326(1) | {1,17,78,1326} | 3 | 3 | **complete** |
| X_0^141(2) | {1,3,94,282} | 3 | 0 | short by 3 |
| X_0^141(2) | {1,6,94,141} | 3 | 1 | short by 2 |
| X_0^77(10) | full | 3 | 0 | short by 3 |

So some of the open bases are directly attackable as the method stands and some are not, and
which is which is a two-second computation. Extending to the oldform part means handling an
eigenvector of `U_p` inside an oldclass, which Theorem 3.1 does not cover.

**The whole open set has now been measured** (`nold_open61.txt`, 2026-10-16): over the 61
quotients the pipeline rerun of PR #117 leaves open, the newforms of level exactly `D N` span the
space at **25** and fall short at 36. Restricting to the 34 that are still open once the three
arguments recorded in `../hyperelliptic/open-cases/` are applied — the metric dual graph of
PR #116 (22), a non-hyperelliptic special-fibre component at a prime exactly dividing `N` (3), and
the composite-level models of `X_0^6(35)` (2, conditionally) — the method reaches **17 of the 34**.
That is a good deal more than the three rows above, which were the only ones measured when they
were written. At `N = 1` the count is automatically complete, and the check is not empty there:
it reproduces the genus at all nine of the `N = 1` open quotients, which is what tests the
Atkin–Lehner convention.

## Stage 3: Noether's count of canonical quadrics

`stage3.m` builds the differentials of `X_0^D(1)/W` from Theorem 3.1 — the weight-2 eigenforms
fixed by `W`, evaluated at `m` points in one enumeration per point — and counts the quadrics
through the canonical image. At genus 3 the count is 1 if the curve is hyperelliptic and 0 if it
is not, so the test is the rank of an `m x 6` matrix of quadratic monomials in the values.

The number of `W`-fixed eigenforms is checked against the genus before anything is evaluated,
which is what verifies the Atkin–Lehner convention rather than assuming it: on the quaternion
side the eigenvalue at `p | D` is minus the classical one, so the eigenvalue of `w_m` is the
product of the classical `a_f(p)` over `p | m`. It came out right at `X_0^91(1)/w_7`, where the
convention does something, as well as at the two `W = {1}` bases where it does not.

| base | genus | forms / genus | points | cutoff | lattice points | seconds | relations | gap | verdict |
|---|---|---|---|---|---|---|---|---|---|
| X_0^35(1) | 3 | 3/3 | 18 | 46 | 1 185 802 | 89 | **1** | 8.1e14 | HYPERELLIPTIC |
| X_0^57(1) | 3 | 3/3 | 18 | 46 | 2 026 839 | 149 | **1** | 1.2e11 | HYPERELLIPTIC |
| X_0^91(1)/w_7 | 3 | 3/3 | 18 | 40 | 1 705 938 | 114 | **0** | — (nothing near zero) | NOT hyperelliptic |

**The ground truth, and a wrong one I started from.** I first took the control pair to be
`D = 35` (hyperelliptic) against `D = 57` (not), from a recollection of Ogg's classification of
hyperelliptic Shimura curves. That recollection was **wrong**, and the repository's own stored
models say so: `data/models/models_57_1.m` records the quotient of `X_0^57(1)` by `w_19` as having
genus 0, so there is a degree-2 map to a rational curve and the curve is hyperelliptic. The
computation agreeing with it is therefore a confirmation, not a contradiction — but only the
model file settles it, and the check was worth making before writing the opposite down.

The genuine non-hyperelliptic control is `X_0^91(1)/w_7`, found by `scan4.m`: it has genus 3 and
every Atkin–Lehner quotient one step above it has genus 1, so there is no degree-2 map to a
rational curve, and for these curves the automorphism group is the Atkin–Lehner group. That is an
argument rather than a table lookup, and it is what the 0 above should be read against.

**A verdict is only taken on a clean gap, and the first attempt did not have one.** At 14 points
clustered in a small region of the upper half plane, `X_0^57(1)` gave a smooth ramp of singular
values — 4.2e-11, 1.1e-7, 1.0e-4, 0.038, 0.42, 2.41 — and reading off the largest ratio named one
relation on what was really no resolution at all, the same defect the Voight–Willis work found in
`GapDimRel`. Nearby points give nearly proportional rows, and a matrix of effective rank near 6
always has a conic through it. Spreading the points by the square root of the arithmetic volume
and raising the cutoff turned that into a ratio of 1.2e11 followed by 52, which is a gap. The
script now refuses a verdict below a ratio of 1e4.

⚠ `stage3.m` itself only does `N = 1`. **`stage3n.m` does general `N`**: it keeps only the
newforms of level exactly `D N`, takes the Atkin–Lehner eigenvalue of `w_m` to be the product of
`a_f(p)` at `p | D` and of `-a_f(p)` at `p | N`, and refuses to run unless their number equals the
genus, so a base the theorem does not reach stops rather than returning a count taken over a
subspace. It reproduces all three controls through `stage3.m`'s own settings to every printed
digit.

### Verdicts on the open quotients (2026-10-16)

Eight of the 34 still-open quotients are decided, every one of them at two or more cutoffs, and
every one **hyperelliptic**: `X_0^39(11)`, `X_0^55(7)`, `X_0^6(185)`, `X_0^46(31)`, `X_0^6(301)`,
`X_0^1326(1)`, and at genus 4 `X_0^2002(1)` and `X_0^2030(1)`. The residuals, gaps and timings are
in `stage3_verdicts.txt`, together with two cross-checks against PR #116 — `X_0^141(1)/w_3` and
`X_0^57(5)` both come out not hyperelliptic, which is what the metric dual graph says — and the
account of the one case that needed more than the cutoff to settle.

That so many come out hyperelliptic is the expected direction, not a surprise: every filter in the
pipeline rules a curve OUT, so what survives to the open set is enriched in curves that really are
hyperelliptic.

### ⚠ The residual floors at the working precision, not at the cutoff

Raising the cutoff stops helping once the truncation error falls below the working precision. At
`X_0^1326(1)` with 30 digits the residual is 9.128e-13 at cutoff 46 and 9.144e-13 at cutoff 50 —
no improvement — while the rest of the spectrum agrees to twelve digits. Raising the precision to
50 digits at cutoff 46 drops it to 8.18e-15, a factor of 112. The sum there is 59 million lattice
terms; at `X_0^35(1)` and `X_0^57(1)`, one or two million terms, the truncation error at cutoff 46
is still the larger of the two and precision changes nothing. **Above roughly 10^7 lattice terms,
raise `PR` along with `CUT`.** This is a limit on the cost law of the previous section, which was
measured at cutoffs 26 to 32 and is right there.

A stalled residual is exactly what a genuine, very small singular value would look like, so it has
to be settled rather than assumed — a residual that moves with the precision is numerically zero,
one that does not is real. Two other explanations were tested first and **refuted**: the fixed
`SCL := 10^16` rounding of the Gram matrix changes no digit of the answer (`nelson3hi.m`,
`stage3hi.m`), because the scaling only decides which vectors the enumeration offers and each is
then accepted or rejected by `4 pi |X| <= cut` at full precision; and widening the
Fourier-coefficient margin `nmax` changes no digit either (`stage3x.m`, `NMX`) — from 20 to 4000 at
`X_0^57(1)`, and from 20 to 20000 at `X_0^1326(1)`, where coefficients to `n = 37769` against the
17789 the cutoff asks for give the same residual to all thirty digits and the same 59 352 084
lattice points. The margin of 20 is ample and nothing is being dropped.

## Files

`wk.m` the Chebyshev machinery for `W_k`; `nelson.m` the first readable implementation, which the
structural checks were run against; `nelson2.m` the working library, with the inner loop reduced
to four complex multiplications and one polynomial evaluation; `nelson3.m` the version for large
bases, enumerating by a `ShortVectorsProcess` and evaluating several eigenforms in one pass.
Drivers: `stage1.m`, `cmpvw.m`, `control.m` (superseded by `control2.m`), `control2.m`, `cost.m`,
`cost2.m`, `stage3.m`, `stage3n.m` (general `N`), `eig.m`, `nold.m`, `scan3.m`, `scan4.m`.
Diagnostics kept because the questions they answer recur: `nelson3hi.m` + `stage3hi.m` scale the
Gram matrix with the precision instead of by a fixed `10^16`, and `stage3x.m` takes the
Fourier-coefficient margin as an argument `NMX`; neither changes any answer, which is the point.
Results: `nold_open61.txt` (how much of each open quotient's space Theorem 3.1 reaches) and
`stage3_verdicts.txt` (the hyperellipticity verdicts, with the residual at each cutoff).

## Traps

* **`W_k` has a singularity at the origin.** A Chebyshev approximation on equal-width panels
  loses twenty digits on the first panel, because the Bernstein ellipse through 0 is what governs
  its convergence there. Panels in a fixed ratio keep the origin equally far from all of them and
  give 1e-38 across the whole range, at 170 times the speed of calling `KBessel` per term.
* **`R^(d)` is not what a Hermite form of the stacked generators returns.** Magma's
  `HermiteForm` on an 8 x 4 matrix gives four rows of the right determinant that span the wrong
  lattice, so the index check `[R^(d):R] = d^2` passes while the lattice is wrong. The symptom is
  a reduced norm with a denominator that cannot occur — 1/4 at `d = 2`, where the only possible
  denominators divide `d`. Use `Lattice(...)` and `BasisMatrix`, and check the denominators, not
  only the index.
* **`QuaternionAlgebra(D)` called twice gives two structures with no covering structure**, so an
  order built against one cannot be coerced into the other. Anything comparing against
  Voight–Willis must pass the same algebra, order and `FuchsianMatrixRepresentation` through.
* **An associative array inside a Magma record is shared with a copy of that record**, so
  spoiling one spoils the other. Running four negative controls in one session made two of them
  return the same number for that reason; `control2.m` gives each its own process.
* **Magma has no `BesselK`**; the modified Bessel function of the second kind is `KBessel(nu, x)`.
* **A control can fire without being independent** — see the `R^(d) = R` row above, which is the
  `d = 1` control multiplied by 2.
* **The transformation law is not evidence about truncation.** Use the rank-one property for that.
* **Clustered evaluation points make Noether's count report a relation that is not there.**
  Rows from nearby points are nearly proportional, and a monomial matrix of effective rank near
  the number of columns always has a quadric through it. The symptom is a smooth ramp of singular
  values rather than a gap, and it is the same symptom a broken reduction produced in
  `../voight-willis/`. Spread the points by the square root of the arithmetic volume, and refuse
  a verdict unless the gap reaches a stated size.
* **Do not take a hyperellipticity verdict from memory of the literature.** Ogg's list, as I
  recalled it, said `X_0^57(1)` was not hyperelliptic; the repository's own model file says its
  quotient by `w_19` has genus 0, so it is. Check `data/models/`.

## hyp16: model fitting, second batch (2026-10-18)

**X_0^55(7)/{1,35,55,77}**, u = w_5: `y² + (x⁴+x+1) y = x⁸ + 5x⁷ + 13x⁶ + 21x⁵ + 24x⁴ + 19x³ + x² − 6x − 2`.
Invariants A = 495/2312, B = −11492/3267, C = −385/228488 recognised at 26 digits from cutoff 46 /
80 digits (238 points; conic residual 5.1e-16, gap 1.7e14), the cutoff-40 run gives the same reduced
minimal model with the roots in the other order (the x → 1/x freedom). Checks: the coset
w_5·W = {5, 7, 11, 385} has 8 fixed points at −20 (w_7) and 8 at −1540 (w_11), fields Q(√5) and
Q(√−35) as the model has over x = 0 and x = ∞; the Jacobian of X/⟨W, w_5⟩ is 385a1 and agrees with the
W-fixed +1 newform at 22 primes; 1 + 2 newform orbits; 62 point counts at 31 primes < 150 and zeta
numerators at all good p ≤ 60 with no mismatch. Still numerical plus point counts, not a proof.
Files: `hyp16/fit_55_7_c46.log`, `fin_55_7_c46.log`, `cert_55_7_c46.log` and the `_c40` trio.

**X_0^6(185)/W**, u = w_3: not a candidate yet. At cutoff 40 the conic is clean (1.7e-13) but only
A = −602/54195 is recognised. Cause is conditioning, not precision: ten of the twelve centres have
|x| in 0.4–1.0 and only two beyond 2.9, so in t = x² the top coefficients rest on two points.
`xscan.m` (cutoff 24, 19 s per point) locates |x| growing toward a pole near z ≈ −4.5 + 0.5i; eight
extra centres (`ctr_6_185_x.txt`, |x|/σ from 2.6 to 6.5) are being evaluated at cutoff 46 / 80
digits and combined with the standard pass by `after6185.sh`. The explicit-centre path is validated
on X_0^35(1) (14 + 2 circles give the control invariants unchanged). ⇒ Check the spread of |x| at
the centres before trusting an unrecognised invariant.

**⚠ `QuaternionOrder(A, N)` is not deterministic across processes**: three runs at D = 55, N = 7
gave two distinct Z-bases for the Eichler order (`detprobe.m`, logs `detprobe1-3.log`), while the
algebra's standard form and the Fuchsian matrix representation were identical. Splitting the
evaluation across processes without fixing the order produces chunks that each lie on a canonical
conic (1e-14) but share none (`conicdiag.m`: 0.10 overall). `hypsetup.m` now writes the order's
Z-basis once and `hypeval_part.m` rebuilds it (discriminant and membership checked); `launch.sh`,
`merge.sh`, `finish.sh`, `queue.sh` drive the chunked pipeline. Control: X_0^35(1) in four chunks
reproduces A, B, C = 1939/952200, 34500/76729, 1939/5000 with the same conic residual 5.17e-16 as
the single-process run (`fit_35_c46_chunked.log`, `fit_35_c46all.log`). `certgen.m` reproduces
every 39_11 claim (`ctl_certgen_3911.log`). ⚠ `zeta3911.m` assembles the zeta numerator with the
functional-equation coefficients reversed; the 39_11 zeta claim rests on `zeta3.m`/`log_zeta3.txt`.

Choice of u (`uchoice.log`): 55_7 w_5; 6_185 w_3; 1326_1 w_2, w_3, w_6 (three independent
instances); 6_301, 6_323, 6_385, 714_5 w_2. Out of reach of this fit: 46_31, 6_379, 6_649 (W is the
whole Atkin–Lehner group, no involution marks the coordinate) and the genus-4 targets (canonical
image a rational normal quartic; y² = degree-10). Per-point cost at cutoff 46 / 80 digits under
load: 55_7 85 s, 1326_1 440 s, 6_301 640 s.

## Oldforms (2026-10-18): `stage3o.m`

The seventeen open quotients that Theorem 3.1 of Nelson does not cover have differentials coming
from newforms of lower level DM, M | N. `stage3o.m` evaluates them: for g new of level DM the
span of g|w_S, S ⊆ T = primes of N/M, with f|w = nrd(w) j(w,z)^{-2} f(wz); w_q for q | DM acts by
the quaternionic eigenvalue, w_p for p ∈ T permutes S ↦ S △ {p} with no scalar (w_p²/p is a unit of
O_N, checked exactly), so the W-fixed differentials are the W-fixed newforms of level DN plus one
v_χ = Σ_S χ(S) g|w_S per (g, χ) with the right product of signs on W. The evaluation is the
level-M lattice sum at the moved point times a phase with the same constant in front, so the
relative scale is fixed by the slash operator. Refuses p² | N, gcd(D,N) > 1, non-squarefree D.
Derivation: `oldforms_derivation.txt`. Validation: the count equals the genus at all seventeen and
48 control candidates (`noldx.m`); the newform path reproduces `stage3n.m` to thirty digits at
X_0^34(7)/w_17 with the same lattice count once both use the same Eichler order (`SEED`); a direct
W-invariance check of every column, with the rejected characters carried along, measures +1 on the
fixed columns and −1 on the rejected ones to 7e-14; twelve stored-model controls across both signs
agree (seven hyperelliptic with residuals 1e-12..1e-14 — including 26_15/{1,3,10,13,30,39,130,390}, all three differentials old across three levels, the two-prime path — and seven not, including the proven smooth
quartics 46_3/w_3 and 34_3/w_2 and the genus-4 57_2/{1,2} with exactly one relation).
First tentative verdict: **X_0^141(2)/{1,6,94,141} NOT hyperelliptic** at cutoff 40 (smallest
singular value 0.21, two of its three differentials old). The rest run from `runq2.sh` /
`items17.txt`; one line per run in `oldforms-logs/queue.log`. Seven of the seventeen are entirely
old (141_2/W1, 57_10, 15_46, 21_34, 77_10, 21_58, 39_34). ⚠ `QuaternionOrder(A, N)` depends on the
random state; use `SEED` when comparing runs.

## hyp16: three defects found while writing the package tests (2026-10-19; package on main's side, PR #127)

* **Galois conjugates are not labelled consistently across processes** (the generator of the
  coefficient field changes), so chunked evaluations of a base with a quadratic newform orbit can have
  the two conjugate columns swapped between chunks. The 6_323 cutoff-46 data had "no conic" for this
  reason while each chunk had one at 1e-17; `fixconj.m` swaps columns 2, 3 in chunks 2 and 3 and the
  conic returns to 1.8e-16 (`ev_6_323_c46f.m`). `hypeval_part.m` now sorts conjugates by their
  coefficients. ⚠ The 6_385 and 714_5 runs were launched before the fix: repair their merged dumps
  with `fixconj.m` before fitting.
* **The base point matters more than expected**: at one Eichler order of X_0^46(3) the default z1 gave
  W-invariance deviations of 1e-3 at cutoff 30 (another base point 1e-9); the package picks z1 among
  six candidates by the smallest relative diagonal value. Mechanism only partly understood.
* **A verdict at a small cutoff can depend on the order**, so tests fix the order by its Z-basis.

**A new case, reported not solved — 6_323.** The repaired values fit perfectly (conic 2.8e-16) but
give B = 22468/134689 rational and A, C complex conjugate (−0.1031 ∓ 0.1124 i): the fixed points of u
on the conic are conjugate over an imaginary quadratic field, the quartic lives there and the Galois
action is x ↔ 1/x. The descent to Q is not implemented. 6_385 and 714_5 may be in the same class.

1326_1 via w_2: invariants −935/324, 1296/3025, −935/324, κ = 17, μ = −34; the model is Q-isomorphic
to the exact one from the height route; check (a) now passes with the per-pair counting of CM fields.
