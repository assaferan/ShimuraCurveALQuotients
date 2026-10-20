# Twisted Weil-polynomial sweep (non-hyperellipticity of X₀^D(N)/W)

Prototype, run locally on 2026-09-25. **The results below are UNREVIEWED.**

## The test

C = X₀^D(N)/W, h an involution of C defined over Q, p ∤ DN. If C is hyperelliptic, so is its
quadratic twist C_h by h (the hyperelliptic involution is central in Aut(C)). Frob_p on H¹(C_h) is
Frob_p∘h on H¹(C), so

    P_{C_h}(t) = P₊(t) · P₋(−t),    P_± = ∏ (t² − a t + p)  over the T_p-eigenvalues a on the
                                         h = ±1 eigenspaces of S₂(DN; W = σ)^{D-new}.

So P_{C_h} must be the Weil polynomial of some genus-g hyperelliptic curve over F_p, that is, a line of
`data/hypg<g>q<p>.txt`. If it is not, C is not hyperelliptic. h = 1 is the pipeline's
own `WeilPolynomial` test and serves as the control.

`tw.m` computes this with modular symbols (sign 0, D-new cuspidal, W-fixed part K). Its asserts are:

* the charpoly of T_p on each eigenspace is a square;
* P_h is monic of degree 2g with constant term p^g;
* the t^{2g−1} coefficient is −Tr(T_p h), the twisted trace of `sweeps/twisted_trace`;
* P_h ≡ P_1 mod 2;
* P_1 is the untwisted polynomial.

## Admissible h

These are the same as for the twisted trace (`sweeps/twisted_trace/README.md`):

* the residual ALs w_Q with Q ∉ W;
* S₂, V₂ and V₃ times an AL, and their pairwise products, under the descent conditions of
  `CheckModularNonALInvolutionModSym`;
* **V₃ only when 9 ∈ W**, because σ(V₃) = V₃W₉. V₃ ops with 9 ∉ W are dropped outright.

An op is kept only if it is an involution on H¹(C) (M² = 1) and commutes with T_p at every tested p.
Otherwise it is listed under `notinv` or `noncomm` and not tested.

## Why parity tests add nothing

P_{C_h}(t) = P₊(t)P₋(−t) ≡ P₊(t)P₋(t) = P_C(t) (mod 2). So any mod-2 or 2-rank criterion gives the
same answer for C_h as for C, and `tw.m` asserts this.

## Table completeness (what the hits rely on)

The user states that the hyperelliptic Weil-polynomial tables are **complete** for:

* g = 3, p ≤ 23;
* g = 4, p = 2, 3, 5;
* g = 5 and 6, p = 2.

Source: the LMFDB (abelian varieties over finite fields, `av_fq_isog.hyp_count`). These are exactly
the (g, p) that `tw.m` and the pipeline use. **Every hit relies on this completeness.** The exhaustive
enumeration below is an independent confirmation, not a prerequisite.

## Coverage

* **Inputs.** `all_in.txt` holds lines of the form `status CurveID D N g DN W`. It has three parts:
  * `controls.txt`: 1355 H and 3 O curves;
  * `known_weil.txt`: 60 curves, `K<p>`, that the pipeline proved by `WeilPolynomial` at p;
  * `undecided.txt`: 584 U curves.
* **Levels.** `levels.txt` has 564 levels and is complete. `levels_ext.txt` has 42 more levels and was
  still running at the snapshot: 5 were in flight and 37 had not started.
* **Primes.** Table primes only: g = 3 at p ≤ 23, g = 4 at p ≤ 5, g = 5 or 6 at p = 2, all with p ∤ DN.
  Curves with g ≥ 7 are not tested.

## Current results (`results_local.out`, the 564 DONE levels)

* **Controls:** 1358 H/O curves gave 0 failures for every h. All 60 K curves fail at h = 1 at the
  pipeline's prime.
* **U curves:** 483 were tested (371 of genus 3, 48 of genus 4, 60 of genus 5, 4 of genus 6). None
  fails at h = 1.
* **Ruled out:** **36 U curves** (`hits.txt`). 32 of them have genus 3 and **4 have genus 4**: 7810,
  9023, 3263 and 3501, all at p = 5. 35 have D > 1, and one has D = 1 (420). These results are
  **UNREVIEWED, and conditional on the completeness of the LMFDB tables** (above).
* **Cross-check:** the 36 failing (curve, w_Q, p) triples with h an AL were recomputed by
  Eichler–Selberg with `TraceDNewALFixed` and no modular symbols (`xcheck.m`, input `xc_all.txt`).
  All 36 match (`xcheck_results.txt`). The V₂ and V₃ hits are not cross-checked.

`hitpolys.txt` lists the distinct failing P_h (20 of them). `vlist.txt` lists the genus-3 ones with
q ≥ 11, taken up to t → −t.

## Completeness check (independent; IN PROGRESS)

This is an exhaustive enumeration of hyperelliptic curves over F_q, run to confirm that each hit
polynomial really is absent. P(−t) is the Weil polynomial of the quadratic twist, so each polynomial
is needed only up to sign.

* **Whole tables** (`tabcheck.m`; `tabcheck_odd.m` in chunks by f(0) for odd q). The table equals the
  enumeration for (g,q) = (3,2), (3,3) and (4,2) (`tabcheck_results.txt`), and for (3,5) as the union
  of the 5 chunks, which gives 1723 = 1723. The run for (4,5) was still going: chunks 0 and 1
  segfaulted in Magma, and chunks 2 to 4 were running.
* **Per polynomial, genus 3, q ≥ 11** (`verify_poly.sh`, which runs `hypsearch.c` and then
  `candcheck.m`). It enumerates the normalised models y² = c(x⁸ + a₆x⁶ + … + a₀) with the same
  (#C(F_q), #C(F_q²)), and Magma computes the exact Weil polynomial of each smooth candidate. None of
  the 10 polynomials in `vlist.txt` was found (`verify_all_snapshot.txt`).
* **Per polynomial, q = 5** (`verify_poly2.sh`, which runs `hypsearch2.c` and `candcheck2.m`, using
  both even and odd models). None of the g = 3 or g = 4 hit polynomials was found. Two rows, those
  with N1 = 3 (g = 3) and N1 = 7 (g = 4), are positive controls that do find curves
  (`verify2_snapshot.txt`).
* Status: **in progress**. `hypsearch.c` and `hypsearch2.c` have not been reviewed either.

## How to rerun

Run from anywhere. Magma must be on the PATH. `.magmarc` is not used (`MAGMA_STARTUP_FILE=/dev/null`).

    sweeps/twisted_weil/run.sh [levels.txt|levels_ext.txt] [JOBS]   # -> out/D_N.out, logs/D_N.log
    python3 sweeps/twisted_weil/summarize.py                        # results_*.out + out/*.out -> hits.txt

A level is complete once its output ends in `DONE`, and complete levels are skipped. To run one level
by hand, from the repo root:

    mkdir -p sweeps/twisted_weil/out
    MAGMA_STARTUP_FILE=/dev/null magma -b D:=6 N:=23 sweeps/twisted_weil/tw.m < /dev/null

The checkers, from the repo root:

    magma -b g:=3 q:=3 sweeps/twisted_weil/tabcheck.m < /dev/null
    magma -b g:=3 q:=5 CHUNK:=0 OUT:=... sweeps/twisted_weil/tabcheck_odd.m < /dev/null
    magma -b IN:=sweeps/twisted_weil/xc_all.txt sweeps/twisted_weil/xcheck.m < /dev/null
    sweeps/twisted_weil/verify_poly.sh 13 '[1,8,59,224,767,1352,2197]'
    sweeps/twisted_weil/verify_poly2.sh 5 4 '[1,6,28,90,230,450,700,750,625]'

The verify scripts compile the C searchers on first use.

Smoke test from this location: (6,23) and (1,30) reproduce the `RES` lines of the local run byte for
byte, and `summarize.py` on `results_local.out` reproduces `hits.txt` exactly.
