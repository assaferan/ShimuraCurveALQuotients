# The even-correction hatch EXISTS at all 28 surveyed obstructed bases — and the cost verdict does not transfer

*2026-09-16. All measurements at `38_5` unless stated. Supersedes the affordability conclusions of
`QUADCONSTRAINTS.md` §8–§10 **only for obstructed bases** — that note's measurements at `34_3` and
`35_1` are correct and unchallenged; what is challenged is extrapolating them.*

## Summary

1. **A two-amount elimination removes the hatch's degree cost — but only at an UNOBSTRUCTED base.**
   Validated exactly at `34_3`; then *proved* not to transfer, and the proof's prediction confirmed
   at `38_5`. Do not re-attempt it on an obstructed base.
2. **The cost figures in `QUADCONSTRAINTS.md` are control-specific.** At `34_3`/`35_1` the price is
   set by which discriminants are integrally solvable. At an obstructed base the amount is *pinned*
   by the charge equation, and at `38_5` the cheapest charge-legal correction costs **2**, not 24.
3. **"Does any even correction exist" is decidable in one shot** by lattice membership, in ~20 s per
   base, replacing a sweep that could only ever see single-discriminant corrections.
4. **The answer is YES at all 28 bases with annihilator data**, using only genuine discriminant
   coordinates. ⇒ conditions 1+2+3 are **jointly satisfiable everywhere surveyed**; condition 3 is
   not a filter on existence. The binding constraints are **condition 4** and **cost**, neither of
   which has ever been measured at an obstructed base.

---

## 1. The two-amount elimination: exact at the control, impossible at an obstructed base

`QUADCONSTRAINTS.md` §6 measures `y²_pert = c_key · f_old · G(s)^amt` with `G` shared across cover
keys. If two amounts `a`, `b=2a` are both legal at one discriminant, then

    (y²_a)² / y²_b  =  (c_a²/c_b) · f_old            -- G cancels EXACTLY

so `f_old` is recoverable at the ORIGINAL degree `2g+2`, from the points the base already supplies,
instead of the `2g+5 + amt·deg Z(disc)` the perturbed `f` demands.

**Measured at `34_3`, disc −164, amounts 6 and 12**, three runs (baseline / 6 / 12), all at the
default `num_vals = 7` (no extra CM points requested):

    branch sets agree across all three runs          42/42 cells
    y²_6²/(y²_12·y²_0) constant per cover key        7/7 keys, EXACT in Q
                                                     (30 usable cells -> 23 independent checks)
    recovered f_old vs baseline after dividing it    exact equality, every cell
    the per-key constant                             a perfect rational square at all 7 keys:
                                                     2/2187, 3/128, 102, 3, 1/51, 1/1003833
                                                     -- each supported on {2,3,17}, the bad primes
                                                     of D*N = 102, exactly as §6 predicts for the
                                                     find_y2_scales row scale
    NEGATIVE CONTROL: y²_6/y²_0 and y²_12/y²_0       take a DISTINCT value at every point
                                                     (5 distinct over 5, etc.) -- the s-dependence
                                                     is real and large; it cancels only in the
                                                     combination

Divisors were checked, not assumed: both perturbed runs print `div_f` = baseline `+ <-164, 6>` resp.
`<-164, 12>`, no mismatch.

### ⇒ WHY IT DOES NOT TRANSFER (this is the part to remember)

At an obstructed base every legal perturbation `P` carries a **fixed nonzero charge**
`φ(P) = -φ(target) ≠ 0`. So:

* two legal perturbations can never be proportional (`P' = λP` with equal charge forces `λ = 1`);
* any multiplicative combination with exponents summing to 1 leaves a residue of charge
  `-φ(target) ≠ 0`, hence a factor of nonzero degree. The residue's cost is bounded below by the
  minimum cost over charge-legal vectors — i.e. exactly the single-perturbation cost.

⇒ **The elimination cannot beat the single-perturbation cost at an obstructed base.** `34_3` allowed
it only because `φ = 0` there (deficit 0): every amount is "in image", which is a control-specific
luxury.

**Prediction recorded before testing, then confirmed at `38_5`**: at most ONE amount per
discriminant is in image. Measured over 12 amounts × 35 discriminants:

    11 in-image hits, at 11 DISTINCT discriminants, exactly one amount each
    disc 4,11,19,500 (phi=2) -> amt 11    disc 20,35 (phi=-1) -> amt -22
    disc 104,264,404,644 (phi=-2) -> -11  disc 520 (phi=-11) -> amt -2
    every disc with 22/phi_j non-integral -> no hit;  phi_j = 0 (disc 24) -> no hit

That also confirms the `φ`-index ↔ `relevant_ds` alignment on 11 independent points, which everything
below depends on.

## 2. The cost verdict is control-specific

`QUADCONSTRAINTS.md` §9–§10 measures the cheapest legal correction as 24 at `34_3` and 48 at `35_1`,
and concludes the hatch is unaffordable by 3×–4.5×. **Both bases are unobstructed controls.** There
the amount is free (`φ = 0`) and the binding filter is integrality, which happens to select expensive
(large `deg Z`) discriminants.

At an obstructed base the amount is not free: `a·φ_j = -φ(target)` pins it. At `38_5`
(`φ(target) = -22`), with `deg Z` measured by `degz.m`:

    disc   phi   amt    degZ   cost = |amt|*degZ   integrally solvable?
      20    -1   -22      1          22            no
      35    -1   -22      1          22            no
     520   -11    -2      1           2            no      <-- cost 2, not 24

⚠ All three fail condition 3, so no single-discriminant correction is usable at `38_5`. That is what
§3 is about. But the *price structure* is completely different from the control's, and **13 of the
35 relevant discriminants at `38_5` have `deg Z = 1`**, so cheap combinations are plentiful:

    degZ at 38_5: 4:1 11:1 19:1 20:1 24:1 35:1 36:1 99:4 100:1 104:3 115:1 120:1 131:5 139:3
    171:4 180:1 196:4 216:3 244:3 264:4 296:5 324:3 340:1 404:7 424:3 456:2 484:3 500:5 520:1
    536:7 596:7 600:8 644:8 676:3      (760 not reached: |d| > 700 sweep bound)

⚠ **`degz.m`'s self-check is hardcoded to `34_3`**, so these numbers are method-validated but not
base-validated, and at an obstructed base there is no committed `div_f` to check the divisor-degree
identity against. Same caveat the `35_1` sweep carries.

## 3. Existence is a lattice-membership question, decidable in one shot

`PROBE_INTSWEEP` tests one discriminant at a time against a guessed amount list, so it can never
answer "does ANY even correction exist" — the charge equation `Σ_j a_j φ_j = -φ(target)` has many
solutions using two or more discriminants and none of those is swept. But writing `L` for the Z-row
span of `dM·coeffs_trunc` and `t` for `dM·target_v`, an even correction exists iff

    t  ∈  L + 2·dM·Z^nds                 (PROBE_MOD2, BorcherdsForms.m)

⚠ **Restrict to the discriminant coordinates.** `Ncols` is consistently `#relevant_ds + 1`; the
trailing column is not a discriminant and the pipeline never perturbs it. Allowing it would let the
test pass using a divisor that does not exist. Both variants are computed; **they agree at all 28
bases**, so the trailing column never carried a witness — but the restricted one is the honest test.

⚠ The test is still a RELAXATION: it permits corrections at ramified discriminants and at CM
evaluation points, which the pipeline avoids for separate practical reasons. **FALSE would be
decisive; TRUE says a witness exists and may not be usable.** It says nothing about condition 4.

### Result: TRUE at 28 of 28

Every base with annihilator data (`annprobe_*.log`): `target_in_L false` (the built-in control — all
are genuinely obstructed) and `target_in_L_plus_2Z_RESTRICTED true`.

    106_3 10_43 10_47 118_3 122_3 142_3 14_17 14_23 14_31 158_3 166_3 22_19 22_23 26_11 26_17
    26_19 34_13 38_11 38_5 46_11 46_7 58_7 62_5 74_7 82_5 86_5 94_3 94_5

### Is the test vacuous?  No — and it is STRICTLY STRONGER than parity at 16 of 28

A check that passes everywhere is not evidence until it is known it could have failed. Two controls:

* **Unit perturbations.** Since `φ` annihilates `L`, membership forces `φ(target)+φ_j` even, so any
  `j` with `φ_j` odd MUST come back false. At `38_5`: 28 true / 8 false, and the 8 are exactly the
  odd-`φ` indices. The test discriminates.
* **Against the parity prediction.** Comparing the measured true-count to
  `#{j : φ^(i)(target)+φ^(i)_j even for all i}`:

      test == parity          12 bases
      test STRICTLY stronger  16 bases   (58_7: 34 measured vs 60 predicted; 74_7: 49 vs 75;
                                          34_13: 30 vs 43; 62_5: 45 vs 53; 82_5/94_3: 37/44 vs 47/54)

⇒ At those 16 bases the TRUE verdict is genuinely new information, **not** a restatement of the
recorded "28/28 `φ(target)` is even" parity survey.

### ⚠ TWO RETRACTIONS OF MY OWN INTERMEDIATE CLAIMS, both the repo's signature failure

* **"22 of 28 bases have 2-torsion in `Z^n/L`" — WRONG, a self-inflicted scaling artifact.** The
  elementary divisors printed were those of `dM·A`, so every base with `dM = 2` showed spurious
  all-2s torsion and `118_3` (`dM = 23`) showed none. Torsion computed that way is meaningless
  across bases. The unit-perturbation control above is scaling-free and is what the conclusion now
  rests on.
* **"The test is just the parity criterion" — WRONG, generalised from one base.** True at `38_5`
  (torsion `Z/9`, odd, so `(Z^n/L) ⊗ F₂` is the free part alone), false at 16 of 28. `38_5` is one
  of the 12 where they coincide.

## 4. What is actually open

    1. EVEN                       satisfiable      (condition 1)
    2. phi(target + P) = 0        satisfiable      (condition 2)
    3. integral solution          satisfiable      (condition 3)  <-- jointly with 1+2, at 28/28
    4. sum_m c(-m) kappa_p(m) in Z at ramified p   UNMEASURED at any obstructed base
    cost = sum_j |a_j| * deg Z(d_j), demand 2g+5+cost rational CM points
                                                   UNMEASURED at any obstructed base

⇒ **Next step: the minimum-cost witness.** Find `a` even with `target + a ∈ L` minimising
`Σ|a_j|·deg Z(d_j)` (a CVP in the weighted norm: take any particular solution from the membership
certificate, then reduce it modulo the lattice of null perturbations `{v even : dM·v ∈ L}`). Then
test that specific correction against condition 4. That produces the first genuine affordability
number for the hatch, and the `38_5` weight table above says the floor is small.

⚠ **One base.** This project's record is that a law fixed from one base gets REPLACED, not refined,
by the second — five times in three days per `QUADCONSTRAINTS.md` §10. Everything in §2 and §4 here
is `38_5` only. The 28/28 existence result in §3 is the only claim here with breadth.

## Reproducing

    # tree: campaign + borcherds-only half of probe-ported-2026-09-12.patch
    #       + ratfit-cmextra.patch + annprobe-obstruction.patch + the PROBE_MOD2 block
    PROBE_MOD2=1 [PROBE_MOD2_CTL=1] magma -b DD:=38 NN:=5 \
        vvdata/weyl-campaign/even-correction/annprobe.m < /dev/null
    PROBE_INTSWEEP=1 PROBE_DS=1 PROBE_SWAMTS=2,-2,4,-4,6,-6,10,-10,11,-11,22,-22 ... (same driver)
    magma -b D_s:=38 N_s:=5 BD:=700 vvdata/weyl-campaign/even-correction/degz.m < /dev/null

⚠ `annprobe.m` does not `SetColumns(0)`, so its output WRAPS MID-TOKEN. Join log lines with `''`,
not `' '`, before parsing, and anchor boolean captures as `(true|false)` — a `\w+` capture silently
swallows the next token and reports `truePROBEANN`. Both bit me.

---

## 5. EXPLICIT WITNESSES AT `38_5` — and the correction is PER COVER KEY

`PROBE_WITNESS` (same patch) enumerates supports of size 1 and 2 with `|amt| <= AMAX` even, testing
membership directly — no reliance on `phi`'s index alignment. At `38_5`, `AMAX = 6`, 15 s per key:

    key 11:  20 witnesses, ALL of support 2   (no singles -- consistent with intsol false on all 11
                                               charge-legal single-disc candidates)
    key 12:  14 witnesses, all of support 2

Costed against `degz.m`'s table (`cost = sum_j |a_j| * deg Z(d_j)`):

    key 11   cost  4   24:+2 520:-2  |  36:-2 180:+2  |  100:-2 340:+2     <-- all have a NEGATIVE amt
             cost 10   20:+2 456:+4                                        <-- cheapest ALL-POSITIVE
             cost 18   104:+4 180:+6  |  120:+6 324:+4  |  120:+6 424:+4
    key 12   cost 18   139:+4 340:+6                                       <-- cheapest ALL-POSITIVE

⚠ **NEGATIVE amounts at a finite discriminant are probably unusable, and that is why the cost-4
witnesses are not the answer.** A negative entry puts POLES in `f`: `f_pert = f_old * G^(-2k)`, so
`f` is no longer a polynomial, while `RationalConstraintsOnEquations` fits
`[s^i : i in [0..2g+2]]`. Extra CM points cannot rescue that — unlike a positive amount, where the
failure is only that the ansatz is too short. (At `34_3` every amount used was positive, and the
baseline's one negative entry sits at disc `-3`, which §7b identifies as the point at infinity.)
**UNTESTED** — if the fit is ever generalised to a rational `f`, the cost-4 witnesses reopen.

⚠ **The correction is PER COVER KEY, not per base.** `coeffs_trunc` is shared but `target_v` is not,
so `phi(target)` differs per key and one vector cannot serve all of them. Applying key 11's witness
makes key 11 pass and the run then stops at key 12 — which `annprobe` had never reached, because it
aborts at the FIRST failing key. ⇒ **the number of obstructed cover keys at a base is not what the
annihilator probe reports; it reports 1.** At `38_5` at least two keys (11, 12) are obstructed, both
have `target_in_L_plus_2Z_RESTRICTED true`, and both have support-2 witnesses.

### What this costs, and what is still unknown

    demand = 2g+5 + cost   rational CM points, per key
    key 11 -> 2g+15        key 12 -> 2g+23

`38_5`'s rational CM supply is **unmeasured and cannot be measured without applying a correction** —
the pipeline dies in `BorcherdsForms`, upstream of `AbsoluteValuesAtCMPoints`, so `#rat` is never
computed at an obstructed base. ⚠ At `34_3` the supply was a hard cap at 10 regardless of the ask;
whether `38_5` behaves that way is unknown. Condition 4 likewise remains untested here.

⇒ **Next: per-key perturbation vectors (`PROBE_EVEN_VEC_<key>`), then a full `genmodels` run at
`38_5`.** That single run answers all three open questions at once: does condition 4 hold, what is
`#rat`, and does the fit close (with `CMEXTRA` if short).

---

## 6. ⚠⚠ THE TARGET COORDINATE IS NOT THE DIVISOR COEFFICIENT — an unguarded hazard in the hatch

**Found 2026-09-16 by making a warning fatal.** Adding `a` to a target coordinate does NOT always
move the divisor by `a`. The ratio depends on the discriminant:

    generic discriminant          factor 1    target +2  ->  divisor +2
    Atkin-Lehner fixed-point      factor 2    target +2  ->  divisor +1
      (|d| = m or 4m, m | D*N)
    d = -3, -4                    factor 4+   target +4  ->  divisor +1
      (extra automorphisms, |Aut| = 6 resp. 4, COMPOUNDING with the AL halving)

Measured at `38_5` (`D*N = 190`): disc 19 (for `w_19`) and disc 20 = 4*5 (for `w_5`) both moved the
divisor by 1 for a requested 2; disc 4 moved it by 1 for a requested 4. The two ramified
discriminants, 4 (`m=1`) and 760 (`m=190`), fit the same rule, which is corroboration rather than
coincidence. Compare [[runaway-class-was-a-scale-bug]], where `d = -4` was the same special case.

### Why this matters more than a wrong cost

**An odd divisor change breaks condition 1** — the cover is no longer preserved, so the "model" is a
DIFFERENT CURVE. At an obstructed base there is no committed model and no oracle, so nothing
downstream would catch it. The whole hatch rests on the correction being even *in the divisor*, and
that was never checked: `probe-ported-2026-09-12.patch` compares `div_f` against a SET union of
`ram` and the requested perturbation, prints `DIVISOR MISMATCH`, and **continues**.

⚠ **The recorded `34_3` results are NOT affected**: disc `-164` is a generic discriminant
(164 = 4*41, 41 does not divide 102), and a re-run this session confirms the change is exactly
`<-164, 6>`. The hazard is latent, not historical — but any future run that selected an AL
fixed-point discriminant would have been silently wrong, and `PROBE_EVEN`'s "prefer the LARGEST
|disc|" heuristic does not avoid them.

### The fix, and why the guard matters more than the rule

`probe-mod2.patch` now (a) excludes `d = -3, -4` outright and requires amounts to be multiples of 4
at AL fixed-point discriminants, with **cost computed from the divisor change** `(amt/factor)*deg Z`
rather than the requested amount; and (b) replaces the print-and-continue check with one that
computes the ACTUAL change `div_f - ram` and **errors** unless every entry is an even integer.

The rule alone would not have been enough: the rule I wrote first covered the AL factor of 2 and
still let disc 4 through at amount 4, and it was the guard — not the rule — that caught it. Keep the
guard even if the rule is later proved complete.

⚠ **Divisor coefficients need not be integers** (half-integers occur at AL fixed points), so the
evenness test must be "is an even integer", not `IsOdd` — which errors on a `FldRatElt`.

### ⇒ CORRECTION to §5: the cheapest cost at `38_5` is 18, not 10

Every cost-10 witness listed in §5 uses disc 19 or 20 at amount 2 — exactly the factor-2
discriminants — so all of them are INVALID: they move the divisor by an odd amount. With those
properly handled the minimum observed cost is **18**, and the fit's demand is `2g+23`, not `2g+15`.
The support-3 search that produced the 10s was costing the requested amount, not the divisor change.
