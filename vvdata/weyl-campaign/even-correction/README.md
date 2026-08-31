# The even-correction escape hatch: measured, and where it is blocked

Work of 2026-08-30/31 on the 28 bases that fail with "Failed to find all Borcherds forms".
Background: [[borcherds-obstruction-is-real]] — the failure is the genuine Borcherds
obstruction, and neither more divisor triples nor deeper poles can ever help. The proposed
escape was: `φ(target)` is **even** and `gcd(φ) = 1`, and a double cover depends on its branch
divisor only mod 2, so an even divisor correction should kill the pairing without changing the
cover.

**Verdict: the precondition holds universally, the correction is constructible, and the hatch
is blocked downstream by an OPEN question in the theory — not by an implementation gap.**

---

## 1. The parity survey — 28/28 EVEN (`annprobe.m`)

The whole hatch turns on one number per base. `annprobe-obstruction.patch` prints the
annihilator `φ` of the divisor map's image and its pairing against the target, at the first
failing key, then aborts — the deficit is invariant across triples (verified at `38_5` over all
96), so triple 1 carries everything and aborting there turns a 900–1700 s exhaustive failure
into a single-key probe. Validated first against the recorded `38_5` result, which it
reproduces exactly (key 11, rank 35/36, ram `[<-760,1>,<-4,-1>]`, `gcd φ = 1`, `φ(target) = -22`)
in **29.8 s instead of 860 s**.

    base    phi(target)  anndim  rank      base    phi(target)  anndim  rank
    38_5    -22          1       35/36     34_13   -28          1       76/77
    10_43   -22          1       76/77     38_11   -2           1       72/73
    10_47   -18          1       84/85     46_11   -20          1       91/92
    106_3   -42          1       64/65     46_7    -22          1       60/61
    118_3   -58          1       71/72     58_7    -46          1       73/74
    122_3   -34          1       74/75     62_5    -48          1       56/57
    142_3   -68          1       86/87     82_5    -10          1       77/78
    14_17   -18          1       43/44     86_5    -24          1       80/81
    14_23   -12          1       32/33     94_3    -40          1       57/58
    14_31   -8           1       81/82     94_5    -32          1       90/91
    158_3   -44          1       98/99     22_23   -22          1       89/90
    26_11   -22          1       50/51     26_17   -36          1       42/43
    26_19   -36          1       85/86
    166_3   -36 / -4     2       99/101    <- 2-dimensional obstruction space
    22_19   -26 / -8     2       41/43     <-
    74_7    -20 / -8     2       92/94     <-

**Every value is EVEN with `gcd(φ) = 1`.** No base is out of reach on parity grounds — the
outcome that would have killed the plan did not happen.

**CORRECTION to [[borcherds-obstruction-is-real]]**, which states the deficit is "exactly 1":
that was measured only at `38_5`. **Three bases have a 2-dimensional obstruction space** —
`166_3`, `22_19`, `74_7`. For those the target must pair to zero against BOTH generators, so a
single shift does not suffice; they need a simultaneous 2-condition solve. Both values are even
in all three cases, so they stay rescuable in principle.

## 2. The positive control at `34_3` (`even-perturbation.patch`)

Obstructed bases have no model to check against, so the premise was tested on a base that
already works and has a committed model. Two design points, both of which mattered:

* **Perturb only the COVER keys, never the hauptmoduls.** `EquationsCovers.m:218` builds
  `d_divs` from the hauptmodul divisors' SUPPORT (keeps `T[1]`, discards multiplicity `T[2]`),
  and `d_divs` drives `CandidateDiscriminants`' `Keep` and `AbsoluteValuesAtCMPoints`' `Include`.
  Perturbing a hauptmodul changes which CM points are selected, so a failure could not be
  attributed to the claim under test — it might be CM starvation.
* **The baseline must reproduce ground truth first.** It does: identical to the committed
  `models_34_3.m` on all 10 keys, plus 5 genus-0 keys the committed file filters out. Worth
  checking, because the pipeline picks an ARBITRARY CM-point pair for the hauptmodul divisors
  (first admitting a rational solution), so different-but-isomorphic presentations were a real
  possibility that would have made "same model" useless as a criterion.

**Result: the even correction IS constructible at the form level.** In every run, all 15 keys
produced a Borcherds form with `div_f` exactly `ram + <disc, 2>`, matching `expected`, zero
mismatches. It then dies downstream in `ValuesAtCMPoints`
(`SchoferFormula.m:1871`, `RationalNumber: s does not represent a rational number`).

## 3. Two refuted explanations, and the diagnostic that settled it

* **REFUTED: collision with the CM evaluation set.** The first run put the correction at
  disc 11, which IS one of the five evaluation points `{3,11,20,24,51}` — the form was made to
  vanish where it is evaluated. Plausible, and wrong: rerunning at discs 180/360, neither an
  evaluation point, failed identically.
* **REFUTED: non-coprimality to N.** 180 and 360 are both divisible by 9 with `N = 3`, which
  fit the known non-rationality axis. Also wrong: discs 164/296, both coprime to 3, failed
  identically. Note "firing" means `N ∤ d_fund`, so discriminants coprime to `N` are exactly the
  FIRING ones — "make it coprime" selects the defective case rather than avoiding it.

`nonrat-diagnostic.patch` then reported the bad cells directly instead of guessing:

    baseline    NONRAT TOTAL  0 cells  (of 9 rows x 7 discs)
    perturbed   NONRAT TOTAL 17 cells, spanning gcd(disc,N) = 1 AND = 3

So the perturbation CREATES the non-rationality (it is not pre-existing), and it is not a
coprimality effect. Bad discs include -11, -56, -68 (coprime) and -51 (not).

## 4. Why it is blocked: an open theory question, not an implementation task

The mechanism is the `KNOWN DEFECT` at `SchoferFormula.m:589`: at a firing discriminant
`Kappa0` returns a log-`N` coefficient of ZERO where it should be `A_m` for `N | m`. A perturbed
divisor acquires components at firing discriminants, whose values are then missing their
`log N` term — so the table entries stop being rational. The `coprime_to_level` filter
(`ShimuraQuotients.m:1420`, self-described as "a blunt instrument" for keeping out "CM points
whose Schofer values misbehave") is the existing workaround for the same defect.

**TESTED AND REFUTED: that `A_m` follows from the preprint's `prop:closedcoef`** (`amtest.m`,
`amtest.out`). That proposition gives the scalar weight-3/2 Eisenstein coefficient `a_E(m)` in
closed form for all `m`, and its hypotheses (`D` squarefree, `ω(D)=2`, `N` prime coprime to `D`)
do cover all 28 obstructed bases and `34_3`. But it reproduces only **1 of the 13** `A_m` values
the code records as measured, and the failure is STRUCTURAL, not a transcription slip:
`prop:closedcoef` carries the embedding support rule, so `a_E(m) = 0` unless every `p | D` is
non-split in `Q(√-m)`. At `15_2, m = 2`: `D₀ = -8`, `χ₃ = (-8/3) = 1`, so `(1-χ_p)` vanishes and
`a_E(2) = 0` — while the measured `A_2 = -1`. Same at `21_2` for `m = 2, 6, 18` (formula 0 vs
measured -2, -4, 2). **`a_E` vanishes exactly where `A_m` does not.**

The reason, from [[b-eisenstein-coefficients-solved]]: the relation is
`A_r = -b^{η*}_0(r)/4`, where `b^{η*}` is the **vector-valued** coefficient at a **nonzero
isotropic coset**, whose support is the level rule `N | r` — a different object from the scalar
`a_E(m)`, whose support is the embedding rule. (The `-b/4` relation itself checks out:
`15_2` has `b(2)=4, b(10)=-4, b(30)=0` → `A = -1, 1, 0`, exactly the code's values.)

**And `b` is an open problem, not merely unimplemented.** That memory records that no product of
local densities reproduces it under any sign convention or assembly variant, even with
Kudla–Yang Props 5.3/5.4 implemented and validated at 470+150 independent checks; and the
follow-on "incoherent series" hypothesis scored 8/9 but was killed by its control.

**Where the preprint stands:** `prop:kappa0` proves the `m = 0` case at nonzero isotropic `μ`
(`κ_μ(0) = -log N/(N-1)`); `prop:closedcoef` proves the all-`m` case for the SCALAR coefficient.
`A_m` needs the intersection neither covers — **general `m` at a nonzero isotropic coset**. That
is the well-posed next theorem, and until it exists the escape hatch cannot be finished.

## Files

    annprobe.m                     the parity probe driver
    annprobe-obstruction.patch     prints phi, phi(target), parity; aborts at first failing key
    annprobe_<base>.log            all 28 survey runs
    even-perturbation.patch        the 34_3 positive control (PROBE_EVEN / _KEYS / _AVOID / _DISC)
    nonrat-diagnostic.patch        reports non-rational cells by row/col/disc/gcd
    control_baseline_diag.log      0 non-rational cells
    control_perturbed_diag.log     17 non-rational cells
    amtest.m / amtest.out          the prop:closedcoef vs measured A_m test (1 of 13)

All three patches are THROWAWAY-branch instrumentation; branches `even-correction` (probe) and
`even-control` (perturbation) are off `619051a`.
