# ppint triage, 2026-08-29 — principal-part integrality across the model backlog

Output of `ppint.m` run over the missing-model backlog on the **production** pipeline
(`tier1-models` at `a924e1a`, i.e. after `origin/main` was merged in).

**`SUMMARY.txt` is the file to read.** It carries the per-base verdicts, the totals, and the
error breakdown. The per-stream `verdicts_*.txt` and per-base `ppint_D_N.log` files are the
raw material behind it.

## What it says

* **Not one previously-untriaged base came back INTEGRAL.** The only INTEGRAL verdicts are
  the three controls (`15_2`, `22_3`, `58_3`), and the three NONINTEGRAL controls
  (`33_2`, `39_2`, `46_3`) also reproduce — 6/6. The earlier belief that integral bases were
  the majority was survivorship: those were bases that had already been worked successfully.
* **`"Failed to find all Borcherds forms"` is the largest failure category**, ahead of
  non-integrality. At `38_5` it is a **rank-1 deficiency** (rank 35 of 36 columns) at a
  single key, identical across all 96 `(infinity, pair)` divisor triples — so the triple
  search at `BorcherdsForms.m:693-700` cannot help. Not a cancellation artifact; unaffected
  by an integral basis or an integral solve.
* CM supply is a **second, independent** blocker: `58_3` has INTEGRAL forms and still cannot
  be generated. Two distinct sites, `BorcherdsForms.m:669` (`< 3` rational CM points — the
  cheapest precondition in the pipeline) and `SchoferFormula.m:1116` (the `max(2g+5)` demand).
  `../cmsupply.m` predicts the second.

## Caveats

* `82_5`'s TIMEOUT is **spurious** — killed as collateral of an experiment. Re-run it.
* TIMEOUTs are at a 2400 s cap and mean *unmeasured*, not *failing*.
* The sweep covered the cheapest bases first; 122 of the 182 untriaged never started.
  `todo.txt` is the full list, `runppint.sh` the driver (run it from the MAIN worktree).
* Every base here was a `polymake_solution_*` **cache hit** (cache covers `M <= 1212`), so
  these results do not depend on the Normaliz backend. It bites only above that frontier.

Untracked by default; commit if wanted. Fuller context in the agent memory note
`handoff-2026-08-29-late`.
