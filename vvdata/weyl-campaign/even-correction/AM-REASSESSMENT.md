# Is the even-correction hatch really blocked on `A_m`? A reassessment

**Date: 2026-09-10.** Reads against `README.md` in this directory (2026-08-31) and
`paper/level-prime-kappa.tex` (2026-09-05, which POSTDATES that README).

> **Verdict: "blocked on `A_m`" is not supported by the README's own evidence.** `A_m` is a
> gauge-fixed repackaging of a correction that the paper now PROVES (`prop:kappa0`) and that the
> code ALREADY IMPLEMENTS (the outer `m=0` term). Implementing `A_m` inside `Kappa0` on top of
> that would DOUBLE-COUNT. The 17 non-rational cells are real and remain unexplained; what is
> withdrawn is the attribution, not the measurement.

This matters beyond the hatch: "the `A_m` theorem" is the stated MAIN LINE in `PLAN.md`, on the
grounds that it unblocks 49 obstructed bases and is the only item that does. That justification
rests on this attribution.

---

## The four steps

### 1. `Sum_m c(-m) A_m = mult(f)` is TRUE BY CONSTRUCTION, not a conjecture about `A_m`

`SchoferFormula.m:589` defines `A_m` as the coefficient of `log N` in `Kappa0(m,d)`, and
[[b-eisenstein-coefficients-solved]] records how the values were obtained:

        mult_f = (1/2) c_eta*(0) = -(1/4) sum_{eta, r>0} c_eta(-r) b^{eta*}_eta(r),
        A_r    = -b^{eta*}_0(r)/4          (the oo-block of that functional)

The `A_m` were SOLVED from that identity (tool `vv/bfit.m`, dimension-0 solution spaces with 5-6
spare conditions). So the identity is the DEFINITION of the values, not a property to be proved
of them.

**Consequence.** Inserting `A_m` into `Kappa0` contributes, at each firing CM point,
`sum_m c(-m) * A_m * log N`, and the per-`m` breakdown is invisible in the CM value, which is a
single number.

⚠ **Be precise about which block that is.** `A_r` is the `oo`-block of the functional; the
`0`-block coefficients `B_j` are a separate family, and `Kappa0` is the `gamma = 0` coset, so an
`A_m`-in-`Kappa0` fix addresses the `oo` block ALONE. The two together are what sum to `mult(f)`:

        mult(f)  =  sum_r c_0(-r) A_r   +   sum_j (0-block) B_j

So the exact statement is: a COMPLETE `A_m`-style fix (`A_r` in `Kappa0` plus `B_j` in
`Kappa(gamma, m)` at `gamma != 0`) contributes exactly `mult(f) log N` -- which is what the outer
term already adds in one step. On the panels where the `oo` block suffices (`rem:gauge`: "the `oo`
block suffices on these panels") the `A_m` half alone already exhausts it. Either way the outer
term is not less complete than the `A_m` route; it is more.

### 2. The code already adds exactly that

`SchoferFormula.m:1024`:

        log_coeffs[i] +:= PointDegree * eta`m0mult * kzero_N;
        //  m0mult  = (1/2) c_eta(0) = mult(f), exact, via M0MultiplierExact
        //  kzero_N = sum_{p | N/(N,d_fund)} log p

This is `prop:kappa0`'s conclusion verbatim: *"the dropped m=0 terms contribute
(1/2) c_eta(0) log N per CM point, i.e. exactly mult(f) sum_{p | N/(N,d)} log p."*

⚠ Note the direction of the dependency. `m0mult` is computed DIRECTLY from the form by
`M0MultiplierExact` — not through any representative of the functional — so it is correct even
for a form whose principal part lies OUTSIDE the panel that pins `A_m`. A perturbed form is
exactly such a form. So on the perturbed input the outer term is, if anything, MORE trustworthy
than an `A_m`-based route would be.

### 3. That code was ACTIVE when the 17 cells were measured

Checked, not assumed. Both hatch branches are off `619051a` (Aug 30 21:12), and

        git show 619051a:SchoferFormula.m | grep -n "M0MultiplierExact\|m0mult\|kzero_N"

returns the outer term at lines 999-1024. So the firing discriminants in the perturbed run
ALREADY received `PointDegree * mult(f) * log N`. The 17 cells are not the signature of that
term's absence.

### 4. The README's §3 diagnostic contradicts its §4 mechanism

§4 states the mechanism as missing `log N` **at firing discriminants**. §3 records the bad discs
as `-11, -56, -68` (coprime to `N`) **and `-51`**. At this base `N = 3`; verified directly:

        d=-11  d_fund=-11  Nprimes=[3]  firing=true
        d=-56  d_fund=-56  Nprimes=[3]  firing=true
        d=-68  d_fund=-68  Nprimes=[3]  firing=true
        d=-51  d_fund=-51  Nprimes=[ ]  firing=FALSE

**And the split is not a stray cell — it is a third of the evidence.** Counting the recorded
`control_perturbed_diag.log` by discriminant:

        disc -68 : 2 cells      disc -56 : 4 cells      disc -20 : 3 cells      disc -11 : 3 cells
        disc -51 : 3 cells      disc -24 : 2 cells

        FIRING cells: 12        NON-FIRING cells: 5   (of 17)

`-24` is the second non-firing disc (`d_fund = -24`, `3 | 24`), and it was not named in the README
prose at all. So **5 of the 17 cells** sit where the level prime divides the fundamental
discriminant.

At `-51` the level prime divides the fundamental discriminant, so `Nprimes` is empty: the outer
term correctly does not fire, AND the `A_m` defect is outside its own stated scope ("at a FIRING
discriminant ... `A_m` for every firing `d`"). **A cell that is bad where no `log N` is expected
cannot be evidence of a missing `log N`.** The README says as much in its own words —
"spanning gcd(disc,N) = 1 AND = 3" — without drawing the consequence.

### 5. And the paper says the support rule is a gauge

`rem:gauge`: *"the rule 'the weight vanishes unless `N | m`' fixes a representative inside that
ambiguity: it is a gauge, not a property of the correction."* This explains, without any further
theorem, the two facts that motivated calling `A_m` open:

* why `N | m` had to be IMPOSED BY HAND for the `bfit` solve to have a unique solution;
* why `prop:closedcoef`'s `-a_E` — a genuine object with the embedding support — disagrees with
  `A_m` pointwise and "reproduces 1 of 13". Two representatives of one functional, differing by
  an element of the panel's relation ideal. Both reproduce every measured multiplier.

⇒ `A_m` is not a quantity appearing in Schofer's formula. It is the amount one would have to
insert into `kappa_0(m)` to reproduce, through the `m`-sum, a correction that in fact lives
entirely in the `m = 0` term at nonzero isotropic cosets. The paper's own `m>0` analysis agrees
that `kappa_0(m)` carries no `log N`: at `m>0` the level prime is never the vanishing place in a
contributing term.

---

## What is therefore WITHDRAWN, and what STANDS

**Withdrawn:** "the next theorem is general `m` at a nonzero isotropic coset, and until it exists
the hatch cannot be finished." The `m=0` case at a nonzero isotropic coset IS the correction, it
is proved (`prop:kappa0`), and it is implemented.

**Stands, unchanged:**
* the parity survey, 28/28 even with `gcd(phi)=1`, and the three 2-dimensional bases;
* the constructibility of the correction at the form level (15/15 keys at `34_3`);
* the measurement `baseline 0 / perturbed 17`;
* the refutation that `A_m` follows from `prop:closedcoef` — TRUE, and now EXPLAINED by
  `rem:gauge` rather than being evidence of a missing theorem;
* `b` being irreducible to a product of local densities of one quadratic space — also expected
  under this reading, since `b` is a gauge-dependent coordinate.

**⚠ A concrete trap this creates.** Implementing `A_m` inside `Kappa0` WITHOUT removing the outer
`m=0` term would apply the correction TWICE. Anyone picking the old plan up cold would do exactly
that, because `SchoferFormula.m:609` says the correction "actually belongs" inside `Kappa0` —
true as bookkeeping, but it is a MOVE, not an ADDITION.

## The better candidate hypothesis

The paper leaves one thing genuinely open in this neighbourhood, and it is narrower. It verifies
the `m>0` local Whittaker VALUES at `N>1` against brute-force representation densities (12 of 12
configurations, including the `p`-scaled odd-`p` shape of `L_-`), and says explicitly:

> "This bounds the reading rather than eliminating it: it verifies the local *values*, not the
> derivative `W'(1)`. But among the `m>0` terms the derivative is only ever taken at a vanishing
> place, and the level prime is never that place in a contributing term."

That last clause is a statement about the divisors the pipeline NORMALLY builds. **A perturbed
divisor is precisely an attempt to build one it normally does not.** If the perturbation creates a
contributing term in which the level prime IS the vanishing place, `W'(1)` is evaluated in a
configuration never validated at `N>1` — and unlike the `log N` story, this predicts failures at
firing AND non-firing discriminants alike, which is what §3 measured.

**This is a hypothesis, not a result.** It is recorded here as the next thing to test, not as a
replacement diagnosis.

## The experiment that decides it

Re-run the `34_3` positive control on CURRENT code with `nonrat-diagnostic.patch`, reporting per
bad cell: disc, firing status, and `eta`m0mult`.

**Predictions recorded BEFORE the run** (per this repo's standing habit):

1. The non-rational cells PERSIST — current code differs from `619051a` here only by the vx fix
   and unrelated work, and the outer term was already present. *(confidence: high)*
2. They INCLUDE non-firing discriminants, reproducing the `-51` class. *(high — it is a re-measurement
   of something already recorded)*
3. At the firing bad cells, `Nprimes` is non-empty and `m0mult` is nonzero, i.e. the outer term
   IS firing and the value is still non-rational. *(medium — if `m0mult` turns out to be 0 on the
   perturbed forms, the picture changes and the outer term is NOT delivering the correction, which
   would partially rehabilitate the old diagnosis in a different form.)*

Outcome 3 is the one that carries information either way, and it is the reason to report `m0mult`
per cell rather than just the disc list.

⚠ The three patches are throwaway instrumentation off `619051a`; they need porting to current
code, and the baseline must be shown to reproduce `models_34_3.m` BEFORE the perturbed run is
believed — that is what made the original control valid and it is not inherited automatically.

---

# RESULTS of the deciding experiment (2026-09-10, same day)

Run in the campaign worktree so the main checkout's running jobs could not compile a half-edited
tree. `nonrat-diagnostic.patch` applied clean; the `even-perturbation.patch` hunk had drifted
(the vx fix and the `IntegralSolution` block moved `BorcherdsForms.m`) and was hand-ported.
Config: `PROBE_EVEN=2 PROBE_EVEN_KEYS=covers PROBE_EVEN_AVOID=3,11,20,24,51`, plus a new
`PROBE_M0` that reports, per CM point and form, whether the outer `m=0` term fired and with what
multiplier.

**The port is faithful**: the perturbed run reproduces the recorded `PROBEEVEN` output exactly --
same perturbation discriminants (296, 164), same `div_f`, same `expected` per key.

## The three predictions, all confirmed

        baseline    NONRAT TOTAL  0 cells   12 populated keys   ok=true
        perturbed   NONRAT TOTAL 18 cells   dies in RationalNumber

⚠ **Not a byte-reproduction of 2026-08-31 (17 cells), and the difference is the CM EVALUATION
SET, not the phenomenon.** That run evaluated at `-56, -68`; this one at `-228, -408`. CM-point
selection shifted between the two code states (the vx fix is the obvious candidate). Same
phenomenon, same scale, different sample -- so quote "17" and "18" as two measurements of one
effect, never as a change in it.

**Prediction 1 (cells persist) -- CONFIRMED.**

**Prediction 2 (they include non-firing discs) -- CONFIRMED, and more strongly than the record:**

        firing      -11, -20                 ->   6 cells
        NON-firing  -24, -51, -228, -408     ->  12 cells

Two THIRDS of the bad cells sit where the level prime divides the fundamental discriminant, and
`PROBE_M0` confirms it from inside the code, per disc: `firing 0 Nprimes [] (no level term is due
here)`. On the 2026-08-31 sample the same split was 12/5.

**Prediction 3 (the outer term fires, with nonzero multiplier, at the firing bad cells) --
CONFIRMED.** At `-11` and `-20` the term fired for all nine forms with nonzero `m0mult`. So at
those six cells the `log N` is PRESENT and the value is non-rational anyway.

⇒ **The `A_m` route cannot repair this.** At 12 of the 18 cells it would add nothing (no level
term is due); at the other 6 it would re-add a term that is already there.

## What the experiment found instead: the perturbed form is NON-INTEGRAL

`m0mult = (1/2) c_eta(0)` at disc `-11`, across the nine forms:

        baseline    -6     3   -3    -9    9   12   -6    6   3      <- integers
        perturbed  -9/4    9  3/4  -21/4  15   18  -9/4   6   3      <- quarter-integers

so `c_eta(0)` goes from even integers to HALF-INTEGERS: the perturbed form is not integral. A
fractional multiplier puts a fractional exponent on a prime (`3^(-9/4)`), which is exactly
`RationalNumber: s does not represent a rational number` -- and, unlike the missing-`log N` story,
it predicts failures at firing AND non-firing discriminants through different primes, which is
what both runs measure.

⚠ **The correlation is strong but NOT exact, and that is recorded rather than smoothed over.** At
`-11` the bad rows are 1, 3, 7 while the quarter-integer forms are 1, 3, 4, 7. Table rows are not
guaranteed to be form indices (the table also carries the `s` and `s~` rows), so the row/form
correspondence must be established before this is called a match.

## The obvious repair, TESTED -- and it fails, which is itself informative

`IntegralSolution` exists for precisely this: `BorcherdsForms.m` notes that `Solution` returns an
ARBITRARY point of `sol + Kernel(coeffs_trunc)`, so an integral representative of the same divisor
may exist and simply not be chosen. If the perturbation merely landed off it, the flag would fix
the hatch outright.

        perturbed + IntegralSolution := true   ->   "Failed to find all Borcherds forms"

So no integral representative of the perturbed divisor exists -- the non-integrality is a property
of the perturbed divisor, not of the solution choice.

**THE CONTROL WAS RUN, AND IT VALIDATES THE VERDICT.** The question was whether
`IntegralSolution := true` is simply too strict at this base, in which case the perturbed failure
would say nothing about the perturbation:

        baseline + IntegralSolution := true   ->   NONRAT TOTAL 0 cells, 12 keys, ok=true

The unperturbed run passes cleanly under the same flag. So the flag is NOT too strict here, and
the perturbed run's "Failed to find all Borcherds forms" is attributable TO THE PERTURBATION.

⇒ **The perturbed divisor admits no integral Borcherds form.** The non-integrality is a property
of the perturbed divisor itself, not of `Solution`'s arbitrary choice within `sol + Kernel`.

⚠ A side observation worth keeping, NOT chased here: the baseline's `m0mult` vector DIFFERS between
the `IntegralSolution` and default runs (`-6 3 -3 -9 9 12 -6 6 3` vs `-6 -3 -3 -9 -9 -12 -6 -3 -3`
at disc `-11`), while both produce 0 non-rational cells and the same 12 keys. So `mult(f)` is not
determined by `div(f)` alone -- the two runs pick forms differing by a trivial-divisor element of
the kernel. Anything that reasons from `mult(f)` as though the divisor fixed it needs checking.

---

# CONCLUSION

**The hatch is not blocked on `A_m`, and it is not blocked on an open theorem.** It is blocked on
INTEGRALITY: the even perturbation, as constructed, produces a divisor with no integral Borcherds
form, whose half-integral `c_eta(0)` puts fractional exponents on primes and destroys rationality
of the CM values at firing and non-firing discriminants alike.

That reframes the hatch from "wait for a theorem" to a SEARCH, and the search is well-posed:
the parity survey already says `phi(target)` is even with `gcd(phi) = 1` at every obstructed base,
so a correction EXISTS in the parity sense. The open question is whether one exists that is also
INTEGRAL.

**Next experiment, cheap and decisive:** sweep the perturbation over (discriminant, even amount)
at `34_3` and test only whether the resulting form is integral -- i.e. whether the run gets past
`BorcherdsForms` under `IntegralSolution := true`. That is a fraction of a full pipeline run per
candidate, and the baseline control above shows the criterion is meaningful at this base. If some
even perturbation IS integral, the hatch reopens on a completely different route from the one
`PLAN.md` currently names.

---

# ⚠ CORRECTION TO THIS FILE'S OWN CONCLUSION (same day, after the sweep)

The CONCLUSION above says the hatch "is blocked on INTEGRALITY". **That is too strong, and the
sweep that was proposed to confirm it is what refuted it.** Integrality is a REAL contributing
cause but NOT the whole cause.

## What the sweep established (and it is real)

Hoisting the test inside the run -- `coeffs_trunc` and `target_v` are already computed, so each
candidate is one linear-algebra test rather than one pipeline run -- swept all 21 discriminants x
{+-2,+-4,+-6} at every key:

    834 candidates:  all in-image (34_3 is unobstructed => surjective divisor map)
                     170 integrally solvable, 664 not

    integral at ALL 7 cover keys, none ramified:
        disc 164  at every amount tested (+-2, +-4, +-6)
        disc  56  at +-4
        disc 180  at +-4

**And it identified exactly what the 2026-08-31 experiment did wrong.** That run's heuristic is
"prefer the LARGEST |disc|", chosen to dodge collision with the CM evaluation set. It therefore
picked disc 296 for keys 8791/8793/8794/8797 and disc 164 for 8792/8795/8796. Measured:

    disc 164, amt +2  ->  intsol TRUE  at every key
    disc 296, amt +2  ->  intsol FALSE at every key it appears

So four of the seven keys were perturbed at a NON-INTEGRAL discriminant, purely because of a
heuristic inside throwaway instrumentation.

## But pinning the integral discriminant does NOT fix the pipeline

Prediction recorded before the run: 0 non-rational cells. **WRONG.**

    PROBE_EVEN_DISC=164, amt +2:   NONRAT TOTAL 11 cells   (was 18)   still dies
    divisors exactly ram + <-164,2> at all 7 keys, no mismatches
    m0mult now -3 9 0 -6 15 18 -3 6 3   <- ALL INTEGERS (was quarter-integers)

So the integrality repair DID work on its own terms -- `c_eta(0)` is integral again and 7 of the 18
cells went away -- and **11 cells remain anyway**, split 4 firing / 7 non-firing. A fractional
multiplier was A cause of non-rationality, not THE cause.

⇒ **The hatch is still blocked, by a residual cause that is NOT identified.** It is not the missing
`log N` (that refutation stands, on its own evidence), and it is not solely non-integrality of the
Borcherds solution. Do not write the next confident single-cause story without a control that could
falsify it -- this file has now produced two.

⚠ What survives unchanged: everything in the "three predictions" section above, the identification
of the 296-vs-164 heuristic error, and the demotion of `A_m` -- which rests on the firing/non-firing
split and on the outer term already supplying `mult(f) log N`, neither of which this correction
touches.

## The residual is INVARIANT — 11 cells, whatever the integral perturbation

| config                     | perturbed?            | NONRAT | m0mult integral |
|----------------------------|-----------------------|--------|-----------------|
| baseline                   | no                    | **0**  | yes             |
| heuristic (296/164, +2)    | yes, 4 keys at 296    | **18** | NO (quarter)    |
| pinned 164, +2             | yes                   | **11** | yes             |
| pinned 164, +4             | yes                   | **11** | yes             |
| pinned 56,  +4             | yes                   | **11** | yes             |
| pinned 180, +4             | NO -- excluded        | **0**  | yes (= baseline)|

⚠ `180` is divisible by `N = 3` and `PROBE_EVEN_COPRIME` defaults to requiring coprimality, so that
run applied NO perturbation at all. Its `m0mult` equals the baseline's exactly. It is therefore an
accidental but useful NULL CONTROL: it confirms the harness reports 0 cells when nothing is
perturbed, so the 11s are caused by the perturbation and not by the instrumentation.

⇒ **Every genuine integral perturbation gives exactly 11 cells, independent of discriminant (164 or
56) and of amount (+2 or +4).** The residual does not depend on WHICH even divisor is added, only
on the fact that one was. That is a structural signal, and it is the single most useful thing this
sweep produced: it rules out "pick a better discriminant" as the remedy, which is precisely what
the previous round of reasoning would have suggested next.

## What to do next, and what NOT to do

⚠ **DO NOT propose a fifth single-cause explanation from this data.** The count so far: CM-set
collision, coprimality to `N`, the missing `log N`/`A_m`, and integrality-alone. Each accounted for
part of the data and was promoted to the whole of it; each was refuted by a control. Integrality is
the only one that survives AS A PARTIAL cause (it accounts for the 18 -> 11 difference, and for the
multiplier becoming integral).

⇒ **The next step is instrumentation, not hypothesis.** Take ONE specific bad cell and print what
`RationalNumber` is actually handed: which prime's exponent is non-rational, and which term of the
Schofer sum contributed it. The bad cells are stable across configurations, so a single cell can be
followed all the way through. That converts guessing into reading, which is what settled the
firing/non-firing question earlier in this file.

---

# THE RESIDUAL, LOCALISED: a fractional exponent on `log 17`, the RAMIFIED prime

`RationalNumber` (LogSum.m:137) fails iff some prime carries a NON-INTEGRAL coefficient -- a
`LogSm` is a formal sum `sum_p coeff_p log p`, so the failure is a FRACTIONAL EXPONENT, not an
irrationality. Naming the prime is therefore the whole diagnostic, and it is cheap.

## The measurement (pinned 164/+2, the 11-cell configuration)

**All 11 bad cells fail on the SAME prime, `p = 17`, with denominator exactly 3.** Coefficients
`+-5/3` and `+-2/3`. EVERY other prime in every bad cell is integral:

    row 3, disc -11  : [<2,-1>,  <11,1>, <17,5/3>]
    row 3, disc -20  : [<2,-2>,  <5,1>,  <17,5/3>]
    row 3, disc -228 : [<3,-2>,  <7,2>,  <17,5/3>, <43,2>]

`D = 34 = 2 * 17`, so **17 is a RAMIFIED prime of `D`** -- not the level prime `N = 3`, where every
previous explanation lived.

## Baseline-vs-perturbed diff at the same cells (63 cells each)

In the BASELINE every `c17` is an integer. The perturbation shifts them, and the shifts are mixed:

    row 2:  -408 -> -2/3        -228/-24/-20/-11 -> +2
    row 3:  -228/-20/-11 -> 11/3        -51 -> +1
    row 5:  -408/-51 -> -2/3           -228/-20/-11 -> +2
    row 6:  -408 -> -8/3
    row 7:  -228/-24/-20/-11 -> 5/3     -51 -> -1

The bad cells are exactly those whose SHIFT has denominator 3; where the shift is an integer the
cell stays fine.

## ⚠ WHAT THIS REFUTES, INCLUDING MY OWN PREVIOUS READING

The natural reading on seeing "fractional coefficient at `p | D`" is that the perturbation breaks
the cancellation which justifies DROPPING the D-part of the m=0 term (`SchoferFormula.m:965`:
"the fractional D-parts cancel against the period / the m>0 Diff-derivatives"). **The arithmetic
does not support it**, and it is recorded as refuted rather than left as a plausible story:

* the D-part at `p = 17` is `(17-1)/(17+1) = 8/9`, which would give denominator **9**; the measured
  denominator is always **3**;
* `1/3` IS `(p-1)/(p+1)` at the OTHER ramified prime `p = 2` -- but that term belongs to `log 2`,
  and `log 2`'s coefficient is integral in every bad cell;
* the m=0 term is per-form times a per-disc INDICATOR, so within a row it would shift all firing
  discs equally. It does not: row 5 shifts `-2/3` at `-408, -51` but `+2` at `-228, -20, -11`.

⇒ The effect is per-`(form, disc)`, not per-form, which points at the **m>0 terms** -- `kappaminus`
at the vanishing place -- rather than at the m=0 constant term. Note `kappaminus` emits a `log p`
coefficient ONLY at a vanishing place, and the perturbation changes the divisor and hence which
places vanish.

## STATUS: localised, not explained

**Established:** the residual is entirely a denominator-3 coefficient on `log 17`, the ramified
prime; every other prime stays integral; the baseline is integral everywhere; the bad cells are
exactly those whose shift has denominator 3.

**NOT established:** the mechanism. Five candidate causes have now been refuted by controls or by
arithmetic (CM-set collision, coprimality to `N`, the missing `log N`/`A_m`, integrality-alone, and
now the dropped D-part). **Do not promote the next pattern to a cause without a control.**

⇒ Next: instrument `kappaminus` itself at one bad `(form, disc)` -- print the vanishing place, the
Whittaker polynomial and `ret` -- since that is the only place a `log 17` coefficient can be
emitted, and the denominator 3 has to enter there.

⚠ Convergence worth noting but NOT leaning on: `PLAN.md` item 0f and
[[b-eisenstein-coefficients-solved]] independently localise the remaining open object to `p | D`.
This arrives at the same primes by a different route. That is corroboration of WHERE, not of WHY.

---

# THE MECHANISM, MEASURED END TO END: an integrality CONGRUENCE on the principal part

Chain, each step measured rather than inferred:

1. `RationalNumber` fails **iff some prime carries a non-integral exponent** (LogSum.m:137). The
   failure is a fractional EXPONENT, not an irrationality.
2. The prime is **always 17**, denominator **always 3** -- and `D = 34 = 2*17`, so it is a RAMIFIED
   prime, not the level prime. Every other prime is integral in every bad cell.
3. It is present BEFORE the final rescaling: `scale = -1/4` (denominator 4, wrong prime), and no
   cell goes integral -> fractional across that step. `c17_pre = -20/3` already.
4. **The principal-part coefficients are INTEGERS: 362 of 362 have denominator 1.**
5. **`Kappa0`'s own log-17 coefficients are FRACTIONAL: 14 denominator 1, 235 denominator 3, 113
   denominator 9.** Ninths are natural at `p = 17`, since `p + 1 = 18` and the ramified factor is
   `(p-1)/(p+1) = 16/18 = 8/9`.

⇒ **The fractional log-17 coefficients are INTRINSIC to `Kappa0` at the ramified prime.** They are
not a defect. In a legitimate Borcherds divisor the weighted sum

        sum_m c(-m) * kappa_17(m)

comes out an INTEGER -- the fractions cancel. The baseline demonstrates it: 0 non-rational cells,
`c17` integral at all 63. **The perturbation breaks that cancellation**, and that is the whole of
the residual.

## What this means for the hatch: there is a FOURTH condition, and it was invisible

A usable perturbation must satisfy, simultaneously:

    1. EVEN                      -- so the double cover is unchanged          (parity survey: 28/28)
    2. phi(target) = 0           -- Borcherds' criterion                      (always solvable, gcd(phi)=1)
    3. integral solution         -- so a Borcherds form exists at all         (measured: 170/834 candidates)
    4. sum_m c(-m) kappa_p(m) in Z  at every ramified p | D and every CM point d   <-- NEW

Condition 4 has been failing silently in every attempt since 2026-08-30. It is not implied by 3:
the `164/+2` run has a fully integral solution AND integral principal-part coefficients, and still
violates 4 at 11 cells.

⚠ **Condition 4 is CHEAP to test**, which is what makes this actionable: `kappa_p(m)` depends only
on `(m, d, Q, lambda)`, not on the model, so a candidate perturbation can be screened by evaluating
one weighted sum per (ramified prime, CM point) -- no covers, no CM tables, no pipeline run. That
is the same hoisting trick that made the integrality sweep cheap.

## Status of this explanation vs the six before it

This is the seventh candidate, so it is worth being explicit about why it is different in kind.
The previous six were HYPOTHESES ABOUT A CAUSE. This is a MEASUREMENT OF WHERE THE FRACTION LIVES,
with counts (362/362 integral coefficients; 348 of 362 fractional kappa values) and with its
control already run -- the baseline satisfies the congruence at all 63 cells, which is exactly what
"0 non-rational cells" means.

⚠ What remains inferred: that condition 4 is SATISFIABLE together with 1-3. Nothing here shows a
perturbation exists meeting all four. That is the next thing to test, and it is now testable
cheaply. If none exists, the hatch is dead for a precise and statable reason rather than an
unexplained one -- which is itself a result.

---

# CONDITION 4 IS SATISFIABLE, AND IT IS A CONGRUENCE ON THE AMOUNT: `amt = 0 mod 2N`

Tested at `34_3` (`N = 3`, so `2N = 6`) with the perturbation pinned at disc 164, which the sweep
showed is integrally solvable at all 7 cover keys for every amount:

    amt    even?  N | amt?   NONRAT cells   fails in
     2      yes     no            11        RationalNumber
     4      yes     no            11        RationalNumber
     8      yes     no            11        RationalNumber
     3      no      yes           23        RationalNumber
     9      no      yes           22        RationalNumber
     6      yes     yes        ** 0 **      QuadraticConstraintsOnEquations
    12      yes     yes        ** 0 **      QuadraticConstraintsOnEquations

⇒ **`amt = 0 mod 2N` clears condition 4 completely** -- 0 non-rational cells, divisors exactly
`ram + <-164, amt>` at all 7 keys, and the pipeline passes `ValuesAtCMPoints` for the first time
since this hatch was attempted on 2026-08-30.

**A PREDICTION WAS RECORDED BEFORE THE RUN AND WAS HALF WRONG, WHICH IS WHY THE LAW IS NOW EXACT.**
Predicted: "condition 4 is divisibility by 3 alone, so `amt = 9` (odd, divisible by 3) will give 0
cells". It gave **22**. Evenness is INDEPENDENTLY required -- odd multiples of `N` are worse than
even non-multiples (22-23 cells vs 11). Had `9` not been tested, "divisible by `N`" would have been
recorded as the law and it is wrong.

**Why the two factors are different in kind**, and the reason the merged statement is the right one:
* EVENNESS is the hatch's founding premise -- a double cover depends on its branch divisor only
  mod 2, so only an even correction leaves the cover alone;
* DIVISIBILITY BY `N` is the new condition -- `kappa_p` at the ramified prime carries denominators
  `3` and `9` (natural at `p = 17`, where `p+1 = 18`), so only a multiple of `N` shifts
  `sum_m c(-m) kappa_p(m)` by an integer.
Together: `amt = 0 mod 2N`.

## The next blocker, which is NOT condition 4

At `amt = 6` and `12` the run now dies later and elsewhere:

    Runtime error in 'QuadraticConstraintsOnEquations':
    Error in Schofer table values at rational points - no solution found!

So the CM values are now rational and usable, and the failure has moved to solving for the
equations. ⚠ **This is a NEW, UNEXAMINED stage** -- do not assume it is the same problem wearing a
different hat, and do not assume it is fatal. It may be CM-supply (the perturbed divisor changes
which rational points are available), which is a known rescue axis with known levers.

## What this implies for the OBSTRUCTED bases -- arithmetic, not yet a claim

At an obstructed base the amount is not free: it must also satisfy `amt * phi(disc) = -phi(target)`.
With `amt = 0 mod 2N` that is a DIVISIBILITY SCREEN, and the parity survey already has every number
needed to run it across all 28 bases without touching a pipeline.

Worked example, `38_5`: `phi(target) = -22`, `N = 5`, so `2N = 10`, and a SINGLE-discriminant
correction needs `amt * phi(disc) = 22` with `10 | amt` -- impossible, since `22` is not divisible
by `10`. ⚠ **That rules out single-discriminant corrections there, NOT the hatch**: a correction
spread over several discriminants has far more freedom, and `gcd(phi) = 1` is exactly the statement
that some integer combination hits any target.

⇒ **Cheapest next thing in the whole file**: run that divisibility screen over the 28 surveyed
bases from the recorded `phi` vectors. No Magma pipeline, no CM tables -- it is arithmetic on data
already in `annprobe_<base>.log`.

---

# THE DIVISIBILITY SCREEN OVER ALL 28 SURVEYED BASES — and a CONFOUND that must be resolved first

If every perturbation component must satisfy `amt = 0 mod M`, then summing
`sum_i amt_i phi(disc_i) = -phi(target)` forces `M | phi(target)`. That is a NECESSARY condition,
computable from the recorded `phi` vectors alone -- no Magma, no pipeline.

## ⚠⚠ THE MODULUS IS CONFOUNDED AT THE ONLY BASE WHERE IT WAS MEASURED

At `34_3`: `N = 3`, and the ramified prime `p = 17` has `p + 1 = 18 = 2 * 3^2`, so the `kappa_17`
denominators (3 and 9) are equally explained by `N` and by the odd part of `p+1`. **The two
readings are indistinguishable there, and I derived "2N" from that single base.**

    reading A:  amt = 0 mod 2N
    reading B:  amt = 0 mod 2*oddpart(p+1)   (p | D ramified)

## The screen under both readings

    reading A (2N):                  3 of 28 pass ->  106_3, 26_11, 82_5
    reading B (2*oddpart(p+1)):      4 of 28 pass ->  10_47, 14_17, 14_23, 62_5

⚠ **The two passing sets are DISJOINT.** So the readings are not a quibble -- they name entirely
different bases, and getting this wrong sends the whole effort at the wrong targets.

## ⇒ The result, stated with its assumption attached

**Under either reading, the hatch as currently understood unblocks 3-4 bases of 28, NOT 49.** That
is a large downgrade of its value, and it is the first quantitative statement of what the hatch can
be worth.

⚠⚠ **THE ASSUMPTION IS UNTESTED AND MAY BE FALSE.** `amt = 0 mod 2N` was measured for a SINGLE
discriminant at a SINGLE base. The real condition is a congruence on
`sum_m c(-m) kappa_p(m)`, which is a statement about the FORM's principal part, not directly about
`phi`. A correction spread over SEVERAL discriminants could satisfy that congruence with individual
components that are NOT separately `0 mod M` -- in which case this screen is far too strict and the
3-4 figure is a floor, not a ceiling. Nothing here tests that.

## The two experiments that settle it, in order

1. **DISAMBIGUATE A vs B.** Re-run the amount sweep at a base where `N` and `oddpart(p+1)` DIFFER.
   `34_7` is ideal: same `D = 34` (so `p = 17`, `oddpart(18) = 9`) but `N = 7`, giving
   `2N = 14` against `2*9 = 18`. It is a working base with a committed model, so it can serve as a
   positive control exactly as `34_3` did. If `amt = 14` clears the cells, reading A; if `amt = 18`
   does, reading B. ⚠ `34_7` was measured at >600 s in the ModelRegen sweep, so budget for it.
2. **TEST THE SINGLE-COMPONENT ASSUMPTION** by perturbing at TWO discriminants at once with
   amounts that are individually not `0 mod M` but sum appropriately. If the cells clear, the
   screen above is void and the hatch's reach is much larger than 3-4.

⚠ Do not quote the 3-4 figure without the assumption. It is a conditional result.

---

# ⚠⚠ `amt = 0 mod 2N` IS REFUTED. The modulus is `D`-dependent, NOT `N`-dependent.

The `34_7` disambiguation was run precisely because `N` and `oddpart(p+1)` coincide at `34_3`.
It killed the law I had recorded.

    base   N   amt   cells   works   2N predicts   mod 6 predicts
    34_3   3    2     11      no        no             no
    34_3   3    4     11      no        no             no
    34_3   3    6      0     YES       YES            YES
    34_3   3    8     11      no        no             no
    34_3   3   12      0     YES       YES            YES
    34_3   3    3     23      no        no             no
    34_3   3    9     22      no        no             no
    34_7   7   28     31      no       YES  <- WRONG   no
    34_7   7   36      0     YES        no  <- WRONG  YES
    34_7   7    4     31      no        no             no

⇒ **`2N` mispredicts BOTH discriminating points at `34_7`**: `28` is divisible by `2N = 14` and
FAILS; `36` is not and CLEARS. **`mod 6` matches all ten measurements across both `N`.**

⇒ **The modulus did not move when `N` went 3 -> 7.** It is a function of `D`, and `N` is eliminated.

⚠ Controls held: `amt = 4` fails at both bases, and all 7 cover keys were perturbed in every run
(`perturbed_keys=7`), so no run is an accidental null like `180/+4` was.

## Still confounded -- three formulas give 6 at `D = 34`

    2*oddpart(2+1)                =  6      (the prime 2)
    2*gcd_{p|D} oddpart(p+1)      =  6      (gcd(3,9))
    2*oddpart(smallest p + 1)     =  6

A base with a DIFFERENT `D` is needed to separate them. Choose one where the three disagree --
e.g. an odd `D`, where there is no `p = 2` to supply the 3.

## ⚠ THE DIVISIBILITY SCREEN IS VOID AND MUST BE REDONE

The screen over 28 bases used `2N`, which is now refuted, so **both of its columns are worthless**:
reading A's `106_3, 26_11, 82_5` and reading B17's `10_47, 14_17, 14_23, 62_5` are NOT candidate
sets. ⚠ **Do not quote either.** The screen has to be re-run once the `D`-formula is pinned, and
until then the hatch's reach is UNKNOWN -- neither the "3-4 of 28" figure nor the original "49"
is supported.

## The methodological point, recorded because it has now happened four times

`A_m`, integrality-alone, `2N`, and before that the `STAR := base_label` reading -- each was fixed
from ONE base and each was underdetermined. In every case the second base did not refine the law,
it REPLACED it. ⇒ **Do not record a law from a single base as settled**, and when two candidate
explanations coincide at the base in hand, say so in the write-up rather than picking the salient
one. `2N` was chosen over `2*oddpart(p+1)` for no better reason than that `N` was the more obvious
number.
