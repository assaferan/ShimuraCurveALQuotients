# `find_t`'s LP: Magma gives up sporadically, and the code read that as infeasibility

Measured 2026-09-15 on Magma **2.29-7** (Mac). `BorcherdsForms.m:find_t` minimises the pole order
of the auxiliary eta-quotient `t` over an integer LP, and bounded every variable below by a
hard-coded **−1000**. It then did

    t, success := Solution(LP);
    assert success eq 0;          // "making sure this is a feasible problem"

so a solver that gives up is reported as `Runtime error in assert: Assertion failed`. That is how
**`6_131` (M = 1572)** and **`6_137` (M = 1644)** were recorded as screen failures.

## They are feasible, with a witness

The three solved neighbours in the same family (`M = 12N`, i.e. `D = 6`) all return the **same**
eta-exponent vector, with only the pole order scaling:

    N =  73   M =  876   eta [0,1,0,-2,-3,6,0,-1,0,2,3,-6]   k = 144 = 2N-2
    N =  89   M = 1068   same eta                            k = 176 = 2N-2
    N = 107   M = 1284   same eta                            k = 212 = 2N-2

Extrapolating `k = 2N-2` to N = 131 and 137 and substituting into `find_t`'s own constraint blocks:

    N 131  M 1572  k 260 : eq true  ge true  le true  ge2 true  ==> FEASIBLE
    N 137  M 1644  k 272 : eq true  ge true  le true  ge2 true  ==> FEASIBLE

⚠ **The −1000 bound was never binding.** The witness's smallest entry is −260. The first hypothesis
— "the bound is too tight at large M" — is **REFUTED by its own witness**, and was drafted and
killed before any edit was applied.

## The failure is sporadic in BOTH directions, so no single bound is safe

    M = 1572   bound -1000                          -> success 25  (gives up)
    M = 1572   bound -261 -300 -500 -800 -1500 -3000 -> success 0, k = 260 EVERY TIME
    M =  732   bound -5000                          -> success 25, while -1000 -1024 -2000 succeed
    M =  948   bound -5000                          -> success 25, while -1000 -1024 -2000 succeed
    M = 1068   bound -5000                          -> success 25, while -1000 -1024 -2000 succeed

Tighter works, looser works, looser-still fails, and the failures do not track `M`. This is a
solver quirk, not a property of the problem — **`success != 0` carries no mathematical
information at all.**

## The fix

`find_t` now retries over `FIND_T_BOUNDS = [-1000, -1024, -2000, -800, -5000, -20000]` and errors
only if every one gives up, saying in the message that this is **not** a proof of infeasibility.

* **−1000 is tried FIRST**, so every base that worked before returns bit-for-bit what it returned
  before. Verified: M = 876/1068/1284 reproduce k = 144/176/212.
* A returned solution whose minimum entry **equals** the bound is rejected and the next bound tried
  — a binding bound may have truncated the search, so that `t` need not be optimal. No case
  observed so far binds; it is a guard, not a code path.

## The test

`tests/FindT.m`, 51 s, 5 bases. It checks the optimum **and** verifies each returned `t` against
`find_t`'s own constraint blocks — because asserting `k = 2N-2` alone would only confirm that the
solver returned what it returned last time, whereas feasibility is checkable independently of the
solver. Negative-controlled both ways, and both controls were **run**, not assumed:

    FIND_T_BOUNDS := [-1000]            -> RED  ("gave up at every lower bound tried", M = 1572)
    ETA perturbed in the last entry     -> RED
    restored                            -> GREEN

## Why this is not a two-base footnote

Of the 78 never-screened even-`D` targets at `#div(M) <= 20`, **24 are `6_N` with M = 12N >= 1788**
— the same family, past the point where the old code first aborted.

⚠ Scope actually measured: `M = 12N` only, Magma 2.29-7. Whether 2.29-10 (lovelace) gives up on the
same pairs is **unchecked**; the retry loop makes that not matter for correctness, but do not quote
the specific failing pairs as version-independent.
