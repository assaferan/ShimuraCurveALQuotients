# The blocker after condition 4 is NOT the quadratic constraints

*2026-09-13. Measured at `34_3`, `amt = 6`, disc 164 — the configuration recorded in
`AM-REASSESSMENT.md` as clearing condition 4 with zero non-rational cells.*

`PLAN.md` and `HANDOFF.md` both record that once condition 4 is satisfied the pipeline clears
`ValuesAtCMPoints` and the blocker **moves to `QuadraticConstraintsOnEquations`**, "a new and
unexamined stage that may be CM supply". This note examines it. The stage attribution is wrong,
the CM-supply guess is right, and the reason it is right is not the one that was guessed.

---

## 1. The error message names the wrong stage — read the `require`, not the traceback

    Runtime error in 'QuadraticConstraintsOnEquations':
    Error in Schofer table values at rational points - no solution found!

That message is raised at `EquationsCovers.m:68`:

```magma
for j->idx in k_idxs do
    B := kernels[j];
    require not IsEmpty(B) : "Error in Schofer table values at rational points - no solution found!";
```

`kernels` is an **argument**. It is built one intrinsic earlier, in
`RationalConstraintsOnEquations` (`EquationsCovers.m:31`):

```magma
M := Matrix([[Rationals()!s^i : i in [0..2*g+2]] cat [Rationals()!rat_y2vals[j]] : j->s in rat_svals]);
B := Basis(Kernel(Transpose(M)));
```

So the failing object is **the kernel of the rational linear fit**, and nothing the quadratic stage
computes is involved — at the moment the error fires, `QuadraticConstraintsOnEquations` has done no
arithmetic at all. `B = 0` says: *the rational CM values admit no polynomial `f` of degree ≤ `2g+2`
with `y² = f(s)`.*

**And the quadratic stage structurally cannot be the blocker.** Its relations are built *from*
`B[1]` and are then solved inside `P(B)` — they select a point in an existing solution space. They
can shrink a positive-dimensional kernel; they can never produce one from an empty kernel.

⚠ At `34_3` the point is doubly moot: **`#quad = 0`.** There are no quadratic CM points at this
base at all, so `relns` is empty and the quadratic stage is a no-op in both the baseline and the
perturbed run. A stage that does nothing cannot be the new blocker.

## 2. What the measurement shows

Two runs of `vvdata/weyl-campaign/genmodels.m` at `D=34, N=3` on current `main` code plus
`probe-ported-2026-09-12.patch`, differing only in `PROBE_EVEN=6 PROBE_EVEN_DISC=164`
(`PROBE_EVEN_COPRIME=0`; **7 `perturb disc` lines**, `div_f` = `expected` at all 7 keys — not a
null run). New instrumentation `PROBE_RATFIT` reports, per cover key, the shape of `M`, its rank,
and the kernel dimension as the degree bound is raised.

    key    W                  g   #ds  #rat  #quad   dimB base   dimB pert
    8791   [1,2,17,34]        1    7     6     0         1           0
    8792   [1,6,17,102]       0    7     6     0         1           0
    8793   [1,3,17,51]        0    7     6     0         1           0
    8794   [1,2,3,6]          1    7     6     0         1           0
    8795   [1,2,51,102]       1    7     6     0         1           0
    8796   [1,3,34,102]       1    7     6     0         1           0
    8797   [1,6,34,51]        1    7     6     0         1           0

**Empty at all seven cover keys**, and the baseline is `dimB = 1` at all seven. At the true degree
bound the perturbed `M` has rank **6 = full row rank** against the baseline's 5: the six values are
genuinely inconsistent with the ansatz, not merely underdetermined by it.

### ⚠ The degree sweep is VACUOUS at this base, and the reason is the whole finding

`#rat = 6` and, for `g = 1`, `ncols = 2g+4 = 6`. The matrix is square. Raising the degree bound adds
columns without adding rows, so rank saturates at 6 and `dimker = ncols - 6` **by counting alone**:

    key 8794 pert:  degbound 4: rank 6, dimker 0     <- the true bound
                    degbound 5: rank 6, dimker 1     <- vacuous
                    degbound 6: rank 6, dimker 2     <- vacuous

⇒ At `MaxNum = 7` the question "does the perturbed data fit a *higher*-degree `f`?" is **not
answerable**, in either direction. Any apparent rescue at a larger degree bound is pure
underdetermination. This is exactly a wrong-object trap: the kernel is nonzero, and it means nothing.

## 3. The hatch's founding premise SURVIVES: the branch locus is unchanged

At every rational CM point, `y² = 0` on the baseline side iff `y² = 0` on the perturbed side:

    branch-locus mismatches: 0 of 42 (rational CM point, cover key) cells

So the even correction does leave the branch divisor alone — the premise the hatch was built on.
What it changes is the `y²` **values at the non-branch points**, and it changes them by a factor
that varies from point to point.

## 4. What the perturbation multiplies `y²` by

Ratios `y²_pert / y²_base` at the non-branch rational CM points factor over a very small set:

    3^10 · 17^8 / 2^10      3 · 17^8 · 43^6 / 2      17^6 / 3^10       31^12 / …
    3^8  · 17^6 / 2^12      17^6 · 43^6 / 24         3^10 · 17^9 / 2^7

The **generic** primes appear only at exponents divisible by `amt = 6` (`43^6`, `31^12`), which is
the signature of the 6th power of the disc-164 form's value. The **bad** primes `2, 3, 17` — that
is, the primes of `D·N = 102`, with `17` the ramified prime of condition 4 — carry the leftover
exponents, and some are **odd** (`17^9`). So the factor is not a square in `Q`, and is not constant.

That is consistent with the structure and needs no new hypothesis: perturbing the divisor by
`amt·Z(164)` multiplies the form by a function with divisor `amt·Z(164)`, whose degree in the base
coordinate is **not zero**. `f` therefore gains degree, and a fit capped at `2g+2` must fail.

## 5. The starvation is a hardcoded constant, not the arithmetic

`#ds = 7` is not what `34_3` can supply. It is

```magma
num_vals := Maximum([2*g+5 : g in genus_list]);   // EquationsCovers.m:337
```

with `max g = 1`, i.e. `2g+5 = 7` — the demand tuned to the **unperturbed** degree `2g+2`
(`2g+4` unknowns, one spare point). `AbsoluteValuesAtCMPoints`' own `MaxNum := 7` default matches.

⇒ **The CM-supply guess in `PLAN.md` is right, but not because the perturbation starves the CM
set.** The CM set is unchanged — same `#ds`, same `#rat`, same six `s`-values in both runs. The
demand is what is wrong: it is computed from a degree the perturbed `f` no longer has.

---

## ⚠⚠ A THIRD NULL-RUN TRAP, same shape as the two on 09-13 — a knob at a dead call site

The first `CMEXTRA` patch went into `EquationsOfCovers` (`EquationsCovers.m:337`, the `[4/6]` path).
Both runs came back with `#ds = 7`, unchanged, in the same 207 s. **`AllEquationsAboveCovers` has
its OWN copy of the `num_vals` computation** (`EquationsCovers.m:~1015`), and that is the one
`genmodels.m` reaches; the `[4/6]` site is not on this driver's path at all. Read naively, the
perturbed run at `CMEXTRA=6` "still failed with more CM points" — a clean refutation of the whole
degree story, **from a run that added no points.**

What caught it: `[4/6]` appears in NO log, baseline included. ⇒ The knob now `printf`s
`CMEXTRA num_vals = %o` unconditionally, so the value it actually used is in the log. Same rule as
`PROBE_EVEN`'s `perturb disc` lines: **print the thing the knob changed, and check it changed.**

## Status of the CMEXTRA experiment

`CMEXTRA=k` (throwaway knob, `EquationsCovers.m`) raises `num_vals` by `k` so the degree question
becomes answerable. Run at `34_3` with `k = 6`, **both perturbed and baseline** — the baseline is
the necessary negative control: if extra CM points break the baseline too, they are bad points and
the perturbed result means nothing.

### The control (read first, and it passes)

    baseline, CMEXTRA=6:  #ds=13  #rat=10  #quad=2   dimB=1 at all 6 keys
    sweep, key 8794 (g=1): deg 4: 10x6  rank 5  dimker 1     <- true bound, now OVERdetermined
                           deg 5: 10x7  rank 6  dimker 1
                           deg 6: 10x8  rank 7  dimker 1
                           deg 7: 10x9  rank 8  dimker 1
                           deg 8: 10x10 rank 9  dimker 1
                           deg 9: 10x11 rank 10 dimker 1
                           deg 10:10x12 rank 10 dimker 2     <- ncols > nrows, vacuous

The extra points are GOOD points: the baseline kernel stays exactly 1-dimensional while the system
is overdetermined, which is what a genuine curve must do (the true `f` of degree 4 also solves every
higher-degree ansatz, and nothing else does). With `nrows = 10` the sweep is non-vacuous through
`degbound = 8`.

### ⇒ PREDICTION, recorded before reading the perturbed log

If the degree story in §4 is right, the perturbed `f` gains `amt · deg Z(164)` in degree, where
`deg Z(164)` is the number of points of disc `-164` in the base's CM cycle — the perturbation adds
`<-164, 6>` to `div_f`, so `f` picks up a factor `(s - s_164)^{6·deg Z(164)}`.

With `amt = 6` and `deg Z(164) >= 1`, the gain is **at least 6**, so for `g = 1`:

    ncols = 2g+4+6 = 12  >  nrows = 10

⇒ **`CMEXTRA=6` is predicted INSUFFICIENT.** Expect `dimB = 0` at every non-vacuous degree bound
(through 8), with the first nonzero kernel appearing only at `degbound >= 10`, where `ncols > nrows`
and it means nothing. **A failure here is therefore NOT evidence against the degree story** — it is
what the degree story predicts. Distinguishing the two readings needs roughly `2g+5+6 = 13`
*rational* points; at the observed 13 ds -> 10 rational ratio that is `num_vals ~ 17`, i.e.
`CMEXTRA ~ 10`.

The reading that WOULD refute the degree story: a kernel appearing at `degbound 6` or `8` (both
non-vacuous). That would say the gain is 2 or 4, not a multiple of 6, and the factor is not
`amt` copies of a CM cycle.

### RESULT: the prediction is confirmed, in both directions

    perturbed, CMEXTRA=6:  #ds=13  #rat=10  #quad=2   dimB=0 at all 6 keys
    sweep, key 8794 (g=1): deg 4: 10x6  rank 6  dimker 0
                           deg 5: 10x7  rank 7  dimker 0
                           deg 6: 10x8  rank 8  dimker 0     <- non-vacuous: gain of 2 REFUTED
                           deg 7: 10x9  rank 9  dimker 0
                           deg 8: 10x10 rank 10 dimker 0     <- non-vacuous: gain of 4 REFUTED
                           deg 9: 10x11 rank 10 dimker 1     <- ncols > nrows, vacuous
                           deg 10:10x12 rank 10 dimker 2     <- vacuous

`M` is at **full COLUMN rank at every testable degree bound**. Degree gains of 2 and 4 — the
readings that would have refuted the degree story — are ruled out on an overdetermined system. The
first kernel appears exactly where counting forces it. As predicted, `CMEXTRA=6` cannot settle the
positive case.

## 6. ⇒ THE DEGREE STORY IS CONFIRMED DIRECTLY, WITHOUT NEEDING A BIGGER RUN

The ratio `h(s) = y²_pert / y²_base`, with the bad primes `{2, 3, 17}` of `D·N = 102` divided out,
is a **perfect 6th power in `Q` at 58 of 58 (point, key) cells** — and its 6th root depends only on
`s`, not on the cover key:

    s        0     1    64/81   4/3   2/3    2     -8     3/4   32/81    1/4
    |G(s)|   1     1     961     1     1     43   4429    713    3013   200777
                        (31^2)                    (43·103) (23·31) (23·131) (41·59·83)

⇒ **`h = c_key · G(s)^6`**, one shared function `G` of the base coordinate, times a per-key constant
supported on the bad primes (that constant is the re-chosen `find_y2_scales` row scale, which is
free to move by a bad-prime factor). `amt = 6` is the exponent. This is exactly what perturbing the
divisor by `6·Z(164)` predicts, and it is measured, not fitted: the exponent 6, the key-independence
of `G`, and the 58/58 exactness were all predicted by the structure before being looked for.

### The consequence that matters for the hatch

`G(s)^6 = (G(s)^3)^2` is a **perfect square**. So

    y² = f_new = c · f_old · (G^3)²      is the SAME CURVE as      y² = c · f_old

under `y ↦ y·G^3`. The even correction does not merely preserve the branch locus (§3) — **it
preserves the whole cover, up to the constant `c`**, which is precisely the quadratic-twist
ambiguity `find_y2_scales` already exists to pin.

⇒ The pipeline is discarding a *correct* answer. It fails not because the perturbed data is wrong
but because `RationalConstraintsOnEquations` solves for `deg f ≤ 2g+2` and the perturbed `f` has
degree `2g+2 + 6·deg G`. Nothing about the mathematics is obstructed here; the ansatz is too short
and the point demand that feeds it is derived from the unperturbed degree.

### What is still unmeasured

`deg G`, hence how far the ansatz must be raised and how many rational CM points that needs:

    need  #rat >= 2g+5 + 6·deg G     (13 at deg G = 1 for g = 1; 19 at deg G = 2)

At the observed `13 ds -> 10 rational` ratio, `deg G = 1` needs `CMEXTRA ~ 10`, `deg G = 2` needs
`CMEXTRA ~ 18`. `deg G` should be `deg Z(164)`, the number of disc-`-164` points in the base's CM
cycle — computable directly, and cheaper than another sweep.

## 7. `deg Z(164) >= 2`, and `34_3`'s rational CM supply is genuinely exhausted at 10

`CMEXTRA=12` (`num_vals = 19`) was run to reach `degbound 10` non-vacuously. It did not:

    #ds = 19   #rat = 10   #quad = 8      (CMEXTRA=6 gave #ds 13, #rat 10, #quad 2)

**`#rat` did not move.** All six extra points are quadratic. So 10 rational CM points are what
`34_3` has — this time a real supply limit, not a hardcoded constant, and it cannot be raised by
asking for more discriminants.

That makes the `g = 0` keys the sharp test, because their ansatz is shorter. For key 8792 the
committed *and* freshly generated `f_old` both have degree **1** (`W=[1,6,17,102]` is the conic
`P![-3,3]`), so `deg f_new = 1 + 6·deg Z(164)`:

    deg Z(164) = 1  =>  deg f_new = 7  =>  ncols = 9 vs nrows = 10    TESTABLE
    observed:  degbound 7: 10x9  rank 9  dimker 0
               degbound 8: 10x10 rank 10 dimker 0

A genuine degree-7 fit would give rank 8, a degree-8 fit rank 9. Both come back at full column rank.
⇒ **`deg Z(164) = 1` is REFUTED, so `deg Z(164) >= 2` and the degree gain is `>= 12`.** Confirming
that directly would need `ncols = 3 + 12 = 15` rational points against the 10 that exist.

⇒ **The direct degree measurement is closed at `34_3` by CM supply.** It does not need reopening:
§6 already establishes the load-bearing claim — the ratio is a perfect 6th power, hence a perfect
square, hence the cover is preserved — without knowing `deg G`.

## 7b. `deg Z(d)` from the divisors — a cross-check that reproduces 7/7 baseline degrees

`f` is a polynomial, so its only pole is at `s = ∞`; the divisor's negative entry is that pole. At
`34_3` every baseline `div_f` has its negative entry at disc `-3`, and

    deg f  =  (multiplicity at -3) · deg Z(3)

    key    W                 div_f (baseline, positive part; pole)      mult   deg f pred   deg f actual
    8792   [1,6,17,102]      <-51,1>;                       <-3,-1>       1         1            1
    8793   [1,3,17,51]       <-408,1> <-24,1>;              <-3,-2>       2         2            2
    8791   [1,2,17,34]       <-408,1> <-51,1> <-24,1>;      <-3,-3>       3         3            3
    8795   [1,2,51,102]      <-68,1> <-24,1>;               <-3,-3>       3         3            3
    8797   [1,6,34,51]       <-408,1> <-68,1>;              <-3,-3>       3         3            3
    8794   [1,2,3,6]         <-408,1> <-68,1> <-51,1>;      <-3,-4>       4         4            4
    8796   [1,3,34,102]      <-68,1> <-51,1> <-24,1>;       <-3,-4>       4         4            4

**7 of 7** with `deg Z(3) = 1` — i.e. the disc `-3` CM point IS the point at infinity. Imposing
`deg(div) = 0` on the same seven divisors then determines the rest, consistently and
overdetermined:

    deg Z(3) = deg Z(24) = deg Z(51) = deg Z(408) = 1        deg Z(68) = 2

(The three independent relations `z_408+z_24 = 2z_3`, `z_408+z_68 = 3z_3`, `z_68+z_24 = 3z_3` are
each confirmed by a second key.) ⚠ Note this is **not** `h(d)/2`: `h(-408) = 4` but `z_408 = 1`,
while `h(-68) = 4` and `z_68 = 2`. Do not fit a formula — read `deg Z(d)` off `FldsOfDefn` as the
sum of the degrees of the fields of definition (`replace_column`, `SchoferFormula.m:1740`, already
uses exactly that quantity), and use this identity as the cross-check.

## 7c. `deg Z(d)` MEASURED at `34_3` — and `deg Z(164) = 4`

`degz.m` computes `deg Z(d)` as the sum of the degrees of the fields of definition
(`FieldsOfDefinitionOfCMPointFast`), the same quantity `replace_column` already uses. It **refuses
to print a sweep unless it first reproduces the five values §7b pins independently** — it does,
5 of 5.

    deg Z(d) at 34_3, |d| <= 500:  43 discriminants
      degZ 1: -3 -11 -20 -24 -27 -51 -75 -147 -228 -267 -312 -408   (12)
      degZ 2: 7    degZ 3: 6    degZ 4: 7    degZ 5: 6    degZ 8: 4    degZ 9: 1

    deg Z(164) = 4

⇒ The correction actually used all along costs a degree gain of **`amt · deg Z(164)` = 6 × 4 = 24**,
so the fit needs `2g+5+24 = 31` rational CM points against a supply of **10**. That is why nothing
downstream of condition 4 has ever worked. (§7's bound `deg Z(164) >= 2`, derived from the rank
refutation alone, is consistent and now superseded by the exact value.)

⚠ **`deg Z(d)` is not `h(d)/2` or any similar formula** — `h(-408) = 4` with `deg Z(408) = 1`,
while `h(-68) = 4` with `deg Z(68) = 2`. Do not fit one; call `degz.m`.

### The cheapest legal correction at `34_3`

Excluding the ramification support (`-408 -68 -51 -24 -3`) and the CM evaluation set
(`3 11 20 24 51`), the `deg Z = 1` candidates are

    -27  -75  -147  -228  -267  -312

At the base's modulus `amt = M = 6` each costs a gain of only **6**, against 24 for disc 164:

    g=1 keys:  ncols = 2g+4+6 = 12  vs 10 rational points   -- still 2 short
    g=0 keys:  ncols =    4+6 = 10  vs 10 rational points   -- EXACTLY determined, TESTABLE

⇒ **The `g = 0` cover keys become testable for the first time.** A perturbation at one of these six
that also satisfies conditions 3 and 4 should give `rank 9, dimker 1` at `degbound 8` on the `g=0`
keys — where disc 164 gives `rank 10, dimker 0`. That is a clean, falsifiable prediction and it
needs no extra CM points.

## 8. ⇒ THE ACTIONABLE CONSEQUENCE: the correction's COST is `amt · deg Z(disc)`

Conditions 1–4 constrain *which* perturbations are legal. None of them mentions the quantity that
decides whether the pipeline can still solve the base afterwards:

    degree gain = amt · deg Z(disc)
    the fit then needs   2g+5 + amt·deg Z(disc)   RATIONAL CM points

This is not a fifth condition — it is a **cost**, and it is the binding one. At `34_3` with
`amt = 6` and `deg Z(164) >= 2` the cost is `>= 12` extra degree, i.e. `>= 19` rational points
against a supply of 10. The hatch is not obstructed there; it is **priced out**.

⇒ **A correction discriminant should be chosen to minimise `amt · deg Z(disc)`.** No selection rule
so far has considered `deg Z(disc)` at all, and `PROBE_EVEN` actively prefers the **largest**
`|disc|` ("least likely to collide with the CM evaluation set") — but large `|disc|` means large
class number and so large `deg Z(disc)`. **That heuristic works directly against the cost.**

The cheapest legal correction is `amt = M` (the base's modulus, 6 at `34_3`) at a discriminant with
`deg Z(disc) = 1`: gain 6, demand `2g+11` rational points. Whether such a discriminant exists among
those satisfying conditions 1–4 is unmeasured — and it is a far better-posed question than the
modulus-formula hunt that was parked.

## 9. MEASURED: at `34_3` the minimum legal cost is 24, against a supply of 10

`costtab.py` crosses `degz.m` against a hoisted `PROBE_INTSWEEP` (`PROBE_SWAMTS=6,12`), over all
21 discriminants in `relevant_ds`. Restricting to those available at **all 9 keys** — a correction
must apply at every cover key or it is not a uniform correction at all — and non-ramified:

    amt = 6                                     amt = 12
    disc  degZ  cost  intsol                    disc  degZ  cost  intsol
     -11    1     6    0/9                       -11    1    12    0/9
     -20    1     6    0/9                       -20    1    12    0/9
     -24    1     6    0/9                       -24    1    12    0/9
     -27    1     6    0/9                       -27    1    12    0/9
     -56    2    12    0/9                       -56    2    24    9/9   <==
     -68    2    12    0/9                       -68    2    24    0/9
    -116    3    18    0/9                      -116    3    36    0/9
    -164    4    24    9/9   <==                -164    4    48    9/9
    -180    4    24    0/9                      -180    4    48    9/9

⇒ **`-164` is the UNIQUE usable discriminant at `amt = 6`** — and it is the joint most expensive
one available. The choice was never free; nothing cheaper is integrally solvable.

⇒⇒ **The minimum cost over ALL legal `(disc, amt)` pairs is 24, reached twice and by different
routes** — `164` at `amt 6` (`6 × 4`) and `56` at `amt 12` (`12 × 2`). Every legal pair measured
here has `amt · deg Z(disc) >= 24`.

    cheapest legal correction at 34_3:  cost 24
    the fit then needs:                 2g+5+24 = 31 rational CM points  (g=1)
    34_3 supplies:                      10

⇒ **The hatch is structurally unaffordable at `34_3`, by a factor of ~3.** Not a bad choice of
discriminant — the cheapest legal one is still 3× over supply. This is a clean negative result
with a stated mechanism, which is worth more than the unexplained failure it replaces.

⚠⚠ **DO NOT PROMOTE "cost >= 24" TO A LAW.** It is ONE base, and this project's record is that a
law fixed from one base gets REPLACED, not refined, by the second (`A_m`, integrality-alone, `2N`,
`STAR := base_label` — four in one day). `24` here coincides with several things at `34_3`
(`4·M`, `2·deg Z(164)·2`, …) and nothing distinguishes them. The honest statement is the
measurement: *at `34_3`, every integrally-solvable `(disc, amt)` pair has `amt · deg Z(disc) >= 24`,
with equality achieved.* Whether a floor exists at all, and whether it is a property of the base,
needs `35_1` or `21_2` — and `degz.m` + `costtab.py` run there unchanged.

### What this does NOT show

* it does not show the hatch is dead in general — only that `34_3`, the positive control, cannot
  afford it. A base with more rational CM points, or a smaller legal cost, is not excluded;
* it does not explain WHY integrality concentrates on large `deg Z(disc)`. The correlation is
  clean here (0 of 11 legal below cost 24) but its mechanism is unexamined, and `deg Z` may be a
  proxy for something else;
* condition 4 was not re-tested for any discriminant other than 164; the sweep measures conditions
  2 and 3 only.

## 10. `35_1`: the mechanism REPLICATES, the cost floor is REFUTED, and condition 4 is DISC-DEPENDENT

### ⚠⚠ RETRACTION: "minimum cost 24" (§9) is a property of `34_3`. `35_1` refutes it.

    base   cheapest INTEGRALLY SOLVABLE disc      cost
    34_3   -164, deg Z = 4                         24     nothing legal below 24 (0 of 11)
    35_1   -112, deg Z = 1, non-ram                12     legal at 2 of 3 cover keys

Fifth time in three days a law fitted from one base has been REPLACED, not refined, by the second.
**What survives is the framework** — `cost = amt · deg Z(disc)`, demand `2g+5+cost` rational CM
points — which came from the 58/58 sixth-power structure, not from the sweep.

### ⚠ AND MY "THIS REOPENS THE HATCH AT 35_1" WAS WRONG, ON TWO COUNTS

* `amt = 6` is **illegal at `35_1`** — the modulus there is 12 — so "cost 6" never existed;
* at the legal `amt = 12`, `-112` **fails condition 4** anyway (dies in `RationalNumber`).

So the cheap discriminant that appeared to break the floor is not usable at all. ⚠ Also: the `g=2`
cover (`9044`, `W=[1,5]`) is **not fitted through this path** — only `9045` (`g=0`) and `9046`
(`g=1`) reach `RationalConstraintsOnEquations` — so a demand quoted "at the `g=2` cover" is a
wrong-object number.

### ✅ THE MECHANISM REPLICATES AT A SECOND BASE

`amt 12` / disc 32 — the recorded clearing config — **reproduces on current code**: condition 4
clears, and the run then dies in exactly the same place as `34_3`:

    RATFIT key 9045 W=[1,35]  g=0  #ds=9 #rat=8 #quad=0  dimB=0
    RATFIT key 9046 W=[1,7]   g=1  #ds=9 #rat=8 #quad=0  dimB=0
    key 9045 sweep: deg 2..6  rank = ncols, dimker 0   (deg 6: 8x8, non-vacuous)
                    deg 7,8   ncols > nrows, vacuous
    Runtime error in 'QuadraticConstraintsOnEquations' ... no solution found

Same error naming the same wrong stage, same empty kernel, same **full column rank**. The §1–§6
analysis is not `34_3`-specific.

### ⇒ CONDITION 4 IS DISC-DEPENDENT — the disc-dependence control, answered

`PLAN.md` asked for exactly this: vary the disc at FIXED base, because base and disc were
confounded (164 at `34_3`, 32 at `35_1`, 16 at `21_2`) and "nothing establishes the modulus is a
property of the BASE rather than of the DISCRIMINANT."

    35_1, amt = 12:   disc  32  ->  condition 4 CLEARS
                      disc 112  ->  condition 4 FAILS

⇒ **Same base, same amount, different discriminant, different outcome.** "The modulus at base X"
is **not a well-defined object**, and every recorded modulus (`M = 6` at `34_3`, `12` at `35_1`,
`| 6` at `21_2`) is really a `(base, disc)` measurement. ⇒ The parked modulus-formula hunt should
be **retired as ill-posed**, not resumed — the six refuted formulas were fitting a quantity that
does not exist.

### Affordability, both bases, using the correction that actually WORKS

    base   working correction     cost      g=1 cols needed   rational supply   short by
    34_3   -164, amt 6           6*4 = 24        30                 10 (capped)    ~3x
    35_1    -32, amt 12         12*4 = 48        54                 12             ~4.5x

⚠ **Supply behaves differently on the two bases**: `34_3` returns 10 whether asked for 7 or 19
(a genuine cap); `35_1` returns 8 when asked for 9 and **12** when asked for 21 (not capped at 8).
So a supply figure is only meaningful with the ask attached.

⇒ **The replacement for the refuted floor is not a number but a COUPLING**: on both bases the
discriminant that satisfies every condition is an expensive one, and the cheap ones are filtered —
by **condition 3** at `34_3`, by **condition 4** at `35_1`. Different condition, same direction.
⚠ Thin: condition 4 has been tested at 2 discriminants at `35_1` and 1 at `34_3`. Do not promote
this to a law either — that is the mistake this section is retracting.

## Summary — what this note establishes

1. **The error message attributes the failure to the wrong intrinsic** (§1). The `require` is on a
   kernel built in `RationalConstraintsOnEquations`. `PLAN.md` and `HANDOFF.md` should be corrected.
2. **The quadratic stage cannot be the blocker** — its relations select inside `P(B)` and can never
   create a kernel — and at `34_3` it has `#quad = 0` and does nothing at all (§1).
3. **The real failure is the rational linear fit**, at full column rank, at all 7 cover keys (§2).
4. **The hatch's founding premise HOLDS**: branch locus preserved, 0/42 (§3).
5. **`h = y²_pert/y²_base = c_key · G(s)^6`** — a perfect 6th power at 58/58 cells once the bad
   primes of `D·N` are removed, with `G` shared across all cover keys (§6). Predicted, then measured.
6. ⇒ **`G^6` is a perfect square, so the COVER is preserved, not just its branch locus** — the
   pipeline is discarding a correct answer because its ansatz is too short (§6).
7. **`deg Z(164) >= 2`**, and `34_3` supplies only 10 rational CM points no matter how many
   discriminants are requested (§7).
8. ⇒ **The binding quantity is a COST, `amt · deg Z(disc)`, not a fifth condition** (§8).
9. **`deg Z(164) = 4`**, so the correction used all along costs 24 and needs 31 rational points
   against a supply of 10 (§7c).
10. **`-164` is the UNIQUE usable discriminant at `amt = 6`**, and the minimum cost over all legal
   `(disc, amt)` pairs at `34_3` is **24** — reached by `164/6` and by `56/12`. Nothing below 24 is
   integrally solvable (0 of 11). ⇒ **the hatch is structurally unaffordable at `34_3`, ~3× over
   supply** (§9). ⚠ One base — do not promote the 24 to a law.

### Reproducing

    export NORMALIZ_BIN=.../normaliz
    git apply vvdata/weyl-campaign/even-correction/probe-ported-2026-09-12.patch
    # plus PROBE_RATFIT in RationalConstraintsOnEquations and CMEXTRA at the
    # AllEquationsAboveCovers num_vals (branch quadconstraints-probe)
    PROBE_RATFIT=1 [CMEXTRA=k] [PROBE_EVEN=6 PROBE_EVEN_DISC=164 PROBE_EVEN_COPRIME=0 \
        PROBE_EVEN_KEYS=covers] magma -b D_s:=34 N_s:=3 OUTDIR=... \
        vvdata/weyl-campaign/genmodels.m < /dev/null
    python3 vvdata/weyl-campaign/even-correction/ratfit_compare.py <baseline.log> <perturbed.log>

⚠ Check `7` `perturb disc` lines and the `CMEXTRA num_vals` line before believing any verdict.

## ⚠ A separate defect noticed in passing: the baseline does not byte-reproduce `models_34_3.m`

The unperturbed `genmodels.m` run gives **15** cover keys against the committed file's 10, and the
shared entries differ. They are **the same curves in a different normalisation**, verified:

    W=[1,2,17,34]   fresh = committed / 9                      (uniform scale, 9 a square: same curve)
    W=[1,6,17,102]  fresh = committed · 9/289                  (uniform scale, a square: same curve)
    W=[1,17]        fresh(x) = committed(17x/3)                (base-coordinate rescale)
    W=[1,102]       fresh(x) = (81/16)² · committed(17x/3)     (both)

So this is the hauptmodul-normalisation axis, not a wrong model — `genmodels.m` does not pin
`base_label` the way the committed generation did. Recorded because a byte-diff against
`data/models/` is **not** a valid reproduction check for this driver, and reading one as a defect
would cost a session.
