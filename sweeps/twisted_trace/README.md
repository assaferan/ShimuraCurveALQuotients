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
primes between 61 and 4g². `report.py` counts the curves affected. For the pending levels
(g = 3..7, so q ≤ 35..195), weil mode computes about 10% more T_p than pb mode (2000 against 1826).

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

## Requirements

* Magma, with `magma` on the PATH.
* This repo checked out at branch `twisted-trace-sweep`.
* GNU `timeout`, from coreutils (on macOS, `gtimeout` is also accepted).
* python3, for the report.

## Run

    git fetch origin && git checkout twisted-trace-sweep      # or clone, then check out the branch
    cd sweeps/twisted_trace
    JOBS=16 TIMEOUT=48h ./run.sh        # restartable: levels whose out/D_N.out has DONE are skipped
    python3 report.py                   # merges results_local.out + out/*.out, dedup by CurveID

Levels run from `levels_pending.txt`, smallest DN first. Each level writes `out/D_N.out`, which is
complete iff it ends in `DONE`, and `logs/D_N.log`. The optional environment variables are `LEVELS=`,
`CURVES=` and `PB=`. One level by hand, from the repo root:

    magma -b D:=15 N:=146 sweeps/twisted_trace/twist.m < /dev/null

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
