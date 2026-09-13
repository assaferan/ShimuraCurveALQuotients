# Probing the `#div >= 24` wall directly, without a pipeline run

*2026-09-13. The backlog plan (`post-m0-backlog-sweep-plan`) excludes 124 missing targets as
"`#div >= 24`, polymake OOM". **polymake is dead** — `nmzsolve.py`/Normaliz replaced it — so that
exclusion rests on a backend that no longer exists and had never been retested. This probes it.*

## The method: call the solver, not the pipeline

`nmzsolve.py M n m out` is the whole polytope solve. Probing it directly costs seconds-to-minutes
per point instead of a multi-hour base build, and it isolates the solver from every other stage.

## ⚠ FOUR TRAPS, all hit or nearly hit while setting this up

### 1. `polymake_script_*` is NOT the cache. `polymake_solution_*` is.

Both live in `polymake/`. The `_script_` files are leftover **polymake input scripts** from the dead
backend (`use application "polytope";`). The real cache is `polymake_solution_<M>_<n>_<m>`, 503
files, read at `BorcherdsForms.m:183`.

I counted `_script_` files and produced a confident table saying **6 of 23** wall-tier targets were
cached. **The true answer is 0 of 23.** Adjacent names, same directory, entirely different object —
and the wrong table would have made the wall test look like it was testing the cache.

### 2. A validated detector for "did Normaliz actually run?"

`nmzsolve.py` always announces itself on **stderr**:

    success:  # M=660 n=1 m=0 k24=12 sq_disc=0 cuspidal=0: 2 lattice points
    timeout:  # normaliz timed out after <NMZ_TIMEOUT>s on <file>.in
    error:    tail of normaliz stdout+stderr, exit 1

A cache hit emits **nothing**. So `grep '^# M='` on a run log distinguishes "solver ran" from
"cache served it" — necessary because `CLAUDE.md` records that above the frontier a fresh solve
fails *silently* ("no solutions", not an error).

**Validated BOTH ways**, which is the point: negative on today's `38_3`/`46_3`/`69_1` builds (0
lines — all three were served entirely from cache and never invoked Normaliz), positive on a forced
fresh solve.

### 3. Reproduce a known value before trusting a new one

First ladder at `M=660`, `n = 1..40`, returned **1-2 lattice points in <1 s each** — which reads as
"no wall" and is meaningless. The cached solutions are large (`276_115_0` is 694 KB), so a real
solve returns many points; 1-2 meant the regime was wrong, not that the problem was easy.

The check that fixes it: re-solve an `(M,n,m)` that already has a cached solution.

    M=276 n=65 m=0:  cached 259 vectors | fresh 259 vectors | IDENTICAL as sets, 1 s

⚠ Compare as **vector sets**, never bytes — committed cache files use escaped `\[` line starts, so
a byte-diff fails on identical content (`CLAUDE.md`).

### 4. `timeout(1)` does not exist on macOS

Use `nmzsolve.py`'s own `NMZ_TIMEOUT` env var (default 1800 s) to bound a probe.

## ✅ RESULT: THE WALL IS REAL ON NORMALIZ. The "do not sweep" routing is CONFIRMED.

Measured at the **matched ratio** `n/M`, which is what makes the two comparable:

    M=276 (#div=12)  n=65   n/M=0.24  ->      1 s   259 points   (reproduces the cache exactly)
    M=660 (#div=24)  n=155  n/M=0.23  ->   >600 s   TIMEOUT
    M=660 (#div=24)  n=100  n/M=0.15  ->    563 s   3 points

⇒ At the SAME `n/M` as a solve that takes one second at `#div=12`, `#div=24` **times out**. The
polytope is nearly empty but ruinously expensive to prove nearly empty. Rungs `n = 250, 400` were
stopped as redundant — they can only be worse.

⇒ **The backlog plan's exclusion of the 124 `#div >= 24` targets STANDS**, on the current backend
and for a different reason than recorded (a time wall, not polymake's OOM). This is a negative
result that saves the sweep rather than enabling it, and it retires a "the wall may be stale"
doubt that would otherwise have kept resurfacing.

⚠ `n = 100` is likely still BELOW the pole order a real `M=660` base needs, so 563 s is a floor on
that rung, not the cost of the base.

⚠ **This probes the SOLVER only.** A base build does many solves plus everything downstream, so
"the solver is slow" is a lower bound on base cost, not an estimate of it.

---

## Corollary (2026-09-13): the wall makes an entire (D,N) CELL unreachable

For odd `D`, `v_2(D) = 0` so `M = 4DN = 2^2 · DN`; with `D` a product of >= 2 odd primes and odd
`N > 1` coprime to it, `#div(DN) >= 8` and therefore

    #div(M) = 3 · #div(DN) >= 24    for EVERY (odd D, odd N > 1)

All 74 such targets are behind this wall by construction. This is why `sweep122` contains only
(even `D`, odd prime `N`) and (odd `D`, `N = 2`) — the missing cell was never a sampling choice.
⇒ The `D`-parity vs `N` confound in the obstructed class cannot be resolved cheaply.
