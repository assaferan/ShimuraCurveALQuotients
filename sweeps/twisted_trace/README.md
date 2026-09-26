# Twisted-trace sweep (non-hyperellipticity of X₀^D(N)/W)

## The test

C = X₀^D(N)/W. Take a prime power q = p^v with p ∤ DN, and an involution h of C defined over Q:
an Atkin–Lehner w_Q with Q ∉ W, or S₂, V₂, V₃ times an AL, subject to the descent conditions of
`CheckModularNonALInvolutionModSym`. V₃ is used only when 9 ∈ W. If C is hyperelliptic, then

    tr := Tr((T_{p^v} − p·T_{p^{v−2}}) ∘ h | S₂(DN; W = σ)^{D-new}) ≥ −(q+1),

so a violation, tr < −(q+1), proves that C is not hyperelliptic. In `twist.m`, the trace is half the
trace on the W-fixed part of the D-new cuspidal modular symbols (sign 0).

**Range of q (`BOUND`).** h is defined over Q, so it commutes with Frobenius. Every eigenvalue of
Frob_q∘h on H¹ therefore has absolute value √q, which gives |tr| ≤ 2g√q. A violation needs
q + 1 < 2g√q, i.e. q < (g + √(g²−1))² < 4g². This is the same bound that `FilterByTrace` uses.

* `BOUND:=weil` is the default. For each curve it tests every q = p^v < 4g², so it is exhaustive:
  no violation can exist beyond that range.
* `BOUND:=pb` reproduces the local run of 2026-09-25: p ≤ 59 and q ≤ 59² for every curve.

These two are **not** nested. The pb range contains more prime powers, but its primes stop at 59, while the
Weil range needs primes up to 4g² (61 at g = 4, and 193 at g = 7). The local results therefore
never went past the Weil bound on q, and all 93 recorded violations satisfy q < 4g² and |tr| ≤ 2g√q
(`report.py` checks this and prints `WEIL-BREACH`). For g ≥ 4, though, the local run skipped
primes between 61 and 4g². The supplement run described below fills that gap. For the pending levels
(g = 3..7, so q ≤ 35..195), weil mode computes about 10% more T_p than pb mode (2000 against 1826).

## The supplement run: why it is needed

The local PB=59 run skipped primes 61 ≤ p < 4g² for curves of genus g ≥ 4. A "no violation" there
was therefore not a complete result. `report.py` computes, for every curve, the set of q that
have actually been tested across all result files. It then compares that set with the Weil range
{q = p^v < 4g², p ∤ DN}. Every U/R curve that has no counted violation and still has a gap goes to
`curves_supplement.txt`. Its level, together with the smallest missing prime as `PMIN`, goes to
`levels_supplement.txt` ("D N PMIN").

At the current snapshot there are **119 curves** (118 U and the reopened curve 9255) at **83 levels**, all
of them small completed levels. Their genera are 4 (47 curves), 5 (67), 6 (4) and 7 (1).

* The missing q are all primes, from 61 to 193. Since 61² > 196 ≥ 4g² for g ≤ 7, no prime
  powers are missing (`report.py` checks this).
* `twist.m PMIN:=61` computes T_p only for those primes. It writes `out/D_N.supp.out`.
* Because the list comes from per-curve coverage and not from a fixed list, it also picks up any
  level that the local queue finishes later in PB=59 mode. To pick those up, copy the queue's newer
  `main.out` in as `results_local.out`; any `results_*.out` is read.
* When everything has run, section 6 of the report ("INCOMPLETE COVERAGE") should read 0.

Coverage counts a q as tested even when a particular op was skipped at q: V₃ ops are only
used at p ≡ 1 mod 3, and an op is skipped if it does not commute with T_p. So "exhaustive" means
exhaustive for the ops the reviewed rules allow.

## Why the numbers are trustworthy

* **Reviewed soundness conditions** (applied in `report.py`):
  * only tr < −(q+1) counts;
  * tr > q+1 is impossible for any curve, so it is flagged as `BUG`;
  * V₃ ops count only when 9 ∈ W;
  * an op that does not commute with T_p at p is skipped;
  * a curve with `noncomm > 0` keeps only its AL/identity violations.
* **Controls:** 914 known-hyperelliptic curves (status H) were tested with 0 violations (913 at review time).
* **Cross-check:** 20 AL-twisted traces agree with Eichler–Selberg (`TraceDNewALFixed`).
* **Independent reproduction:** the reviewer's separate implementation reproduced the violations on
  the 3 reopened curves: 5124 (q=7, −12), 7923 (q=13, −18) and 8387 (q=7, −12).
* **Same code:** `twist.m` is the exact code of the local run, with only its I/O changed. With
  `BOUND:=pb` it reproduces the `RES` lines of `results_local.out` byte for byte at (26,45), and it
  runs without `.magmarc` (`MAGMA_STARTUP_FILE=/dev/null`).
* **Supplement check:** with `PMIN:=7` at (26,45), it finds the known violation on 8387
  (q=7, V3*w1, −12) again. A smoke test at (1,84) and (38,5) with `PMIN:=61` merged correctly and
  dropped the gap count from 119 to 117.

## Requirements

* Magma, with `magma` on the PATH.
* This repo checked out at branch `twisted-trace-sweep`.
* GNU `timeout`, from coreutils (on macOS, `gtimeout` is also accepted).
* python3, for the report.

## Run

    git fetch origin && git checkout twisted-trace-sweep      # or clone, then check out the branch
    cd sweeps/twisted_trace
    JOBS=16 TIMEOUT=48h ./run.sh        # phase 1 pending (Weil), then phase 2 supplement; restartable
    python3 report.py                   # merges results_*.out + out/*.out per CurveID

`run.sh` runs in two phases:

1. The levels in `levels_pending.txt`, smallest DN first. Each writes `out/D_N.out`.
2. The supplement. `report.py` first regenerates the supplement lists, then each level is run with
   its `PMIN` and writes `out/D_N.supp.out`.

A level's output is complete iff it ends in `DONE`, and complete levels are skipped on a rerun. A
partial output is kept as `*.partial-<time>.out`, and its finished `RES` lines still count. Logs go
to `logs/`, and the final report goes to `logs/report_final.txt`.

The phases can be run separately, since the supplement levels are small and quick:

    PHASE=2 JOBS=16 ./run.sh             # supplement only
    PHASE=1 JOBS=16 TIMEOUT=48h ./run.sh # pending only

The other optional variables are `LEVELS=`, `CURVES=` and `SUPPLEVELS=` (subsets, for testing) and
`PB=` (a prime cap). One level by hand, from the repo root:

    magma -b D:=15 N:=146 sweeps/twisted_trace/twist.m < /dev/null
    magma -b D:=38 N:=5 PMIN:=61 IN:=sweeps/twisted_trace/curves_supplement.txt sweeps/twisted_trace/twist.m < /dev/null

To bring the results back, copy `sweeps/twisted_trace/out/` into the same place in the local
checkout and rerun `python3 report.py`.

## Cost and state

The pending set is 264 U/R curves at 146 levels, DN 2190–15330 (`curves_pending.txt`). It was
taken from a `results_local.out` snapshot while the local queue was still running.

Locally, the largest completed levels took 2+ hours each. The first pending levels, (15,146) and
(10,231), had already been running for more than 4 h locally when the snapshot was taken. The cost
is roughly one operator solve per AL divisor plus one per prime (about 2 s per T_p at dim 200, and
growing with the dimension), so expect days at the top end. Use `TIMEOUT` generously.

The current local result is 71 of 563 tested U curves ruled out (70 with D>1, 1 with D=1), with
0 control violations.
