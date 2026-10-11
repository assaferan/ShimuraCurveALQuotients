# The lovelace / lava launch scripts

Copies of the shell scripts that drive model runs on the servers, banked here on 2026-10-11 because
until then they lived only in `~/gymodels/` on lovelace and `/scratch/home/assaferan/` on lava.
`README.txt` is lovelace's per-directory ownership log at the same date (who launched what, where
its logs are, whether it is safe to kill).

| script | what it does |
|---|---|
| `launchM.sh <tree> <logdir> <outdir> <maxp> <M> <n> D:N ...` | runs the bases of one level, at most `maxp` at a time, once `polymake/polymake_solution_<M>_<n>_0` and `polymake/tshift_w0_<M>.txt` are in `<tree>` |
| `launch660.sh` | the same for the nine M = 660 bases, written before `launchM.sh` generalised it |
| `bank_launch.sh <tree> <M> <n>` | enumerates the first m = 0 rung `(M, n, 0)` of a level and its weight-0 shift set, with 10-day and 2-day limits |
| `probe_next.sh` | 30-minute probes, one base per level, to read the first polytope a level asks for off `polymake/nmzsolve.err` |

## Two things the scripts cannot show

**`NMZ_TIMEOUT=7200` in `launchM.sh` is a budget per polytope, and a timeout used to be silent.**
Until PR #81, `BorcherdsForms` read a Normaliz timeout as an empty polytope, so a base whose rung
needs more than two hours climbed the ladder one rung per timeout. At level 924 that cost 3.5 days
per odd-D base before the n-cap fired. With #81 the run stops at the first timeout with the
polytope named; enumerate it separately (`nmzsolve.py M n m out` with a long `NMZ_TIMEOUT`, as
`bank_launch.sh` does), drop the file into `<tree>/polymake/`, and relaunch.

**Odd-D bases need two more polytopes than even-D ones, and the m = 0 shift bank does not cover
them.** For odd D the search also builds the ring of weight-1/2 forms with poles at both cusps,
which at every finished odd-D base converges after exactly two rungs `(M, n0, k)` and
`(M, n0 + k, k)`, k the pole order of the level's `t`: (420, 42, 104) and (420, 146, 104);
(924, 89, 236) and (924, 325, 236). Bank both before launching an odd-D base of a 24-divisor level.

**The queue script outlives the jobs it launched.** `launchM.sh` waits for a free slot and starts
the next base; killing a running base frees a slot and the queue fills it. To stop a level, kill
the `bash launchM.sh` process first, then the Magma jobs (by PID, never `pkill -f magma.exe`).
