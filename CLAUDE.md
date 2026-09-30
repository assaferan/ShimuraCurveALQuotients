# Working in this repo

`README.md` covers loading the curve data. This file covers the things that are true regardless of
what you are working on, and that have each cost a session at least once.

## Start here

> ### ⇒ Spend your effort on WHICH OBJECT a claim is about, not on whether the computation is right.
>
> Learned expensively on 2026-09-04: nearly every wrong conclusion that day had correct arithmetic
> about the wrong object — a rank over monomials when the claim was about forms, a parser emitting
> valid polynomials from truncated input, greps that read a fragment and generalised. Validations
> that check "is this number right" cannot catch "is this the right number". The two habits that
> did catch things: **reproduce a KNOWN value before trusting a new one**, and **draft an edit
> rather than applying it** — a wrong "the paper is wrong" claim was retracted while writing the
> diff. See `HANDOFF.md`, "READ THIS FIRST".

* **`PLAN.md`** — what to do next. Opens with a "Picking this up cold" block.
* **`HANDOFF.md`** — what happened. Authoritative on state; it wins over `PLAN.md` and over any
  agent memory when they disagree.

## Changes reach `main` only through a reviewed PR

`main` changes only when a human (@sachihashimoto or @assaferan) merges a pull request on GitHub.

* **Branch first.** Before the first edit of a task, create a topic branch off an up-to-date
  `main`, in a worktree: `git fetch origin && git worktree add -b <topic> worktrees/<topic>
  origin/main`. If you find changes sitting on `main`, move them to a topic branch before
  committing.
* **One topic per branch**, and commits and pushes go to that branch only. Pushing to `main`,
  force-pushing it, and merging PRs are the human reviewers' job.
* **Open a PR for review** with `gh pr create --base main` when the work is ready, then leave it
  for a human to review and merge. The description states:
  * what changed and why (the mathematical claim or bug being addressed);
  * how it was verified: the exact Magma scripts/tests run and their outcome;
  * every file under `data/` that was regenerated, since test runs can rewrite `data/*.dat` as a
    side effect.
* **Stacked work:** a branch that depends on an unmerged PR branches from that PR's branch and
  targets it with `--base`, so each PR shows only its own diff.
* Merged topic branches are deleted; they are not part of the long-lived set below.

## Running Magma

    AttachSpec("ShimuraQuotients.spec");     // always, before anything else

    magma -b run_tests.m < /dev/null                        # whole suite
    magma -b target:=ModelChecks run_tests.m < /dev/null    # substring filter
    magma -b filename:=tests/Kappa0.m run_tests.m < /dev/null
    magma -b D:=15 N:=2 myscript.m < /dev/null              # args are name:=value

**Always redirect stdin.** On any runtime error Magma drops into its interactive interpreter and
blocks on stdin forever. Piping to `head`/`tail` makes this worse, not better — one such hang ran
4 h 23 m at 0.02 s of CPU. Use `< /dev/null` and redirect output to a file.

**Magma buffers stdout when it is a file.** An unchanged log is *not* evidence that nothing is
happening — one run sat at 154 bytes for hours, then flushed 353 KB at once. A killed run loses
whatever is still in the buffer, so a `kill`-based timeout can leave a 0-byte log that looks
identical to "never started". macOS has no `timeout`; GNU `timeout` (on the Linux boxes) truncates
cleanly.

**Do not edit anything under `tests/` while a suite is running.** `run_tests.m` globs the file list
once at start but *loads each file when it reaches it*, so a run in flight can read a half-written
or deliberately-broken version and report a failure that means nothing. This bit a session on
2026-09-23: negative controls, which rewrite a test file and restore it seconds later, were run
against a suite launched ten minutes earlier, and the whole 30-minute run had to be discarded.
Develop controls in the scratchpad and install them afterwards. It is the same hazard as the
"never `git pull` a clone with jobs running from it" rule below — the working tree is a clone too.

**⚠ THE FULL SUITE DOES NOT COMPLETE ON THIS MAC — it dies at `X0_206_1`.** Measured twice on
2026-09-23, independently: two runs reached 56 and 57 files and both were killed by macOS for memory
pressure at exactly `X0_206_1.m`, with 0 failures up to that point. So a local `run_tests.m` with no
target is a ~3-hour way to learn nothing past the `X0_1*` range. Either use `target:=` /
`filename:=` for the tests you care about, or run the suite on **lava**, which is where the offline
tests already go. ⚠ A killed suite LOOKS like a clean one — check the file count and that
`Tests failed:` is present, per the truncation rule.

**Kill Magma by PID, never `pkill -f magma.exe`.** This Mac hosts several Claude sessions and the
`core` repo sessions run their own `magma.exe`; the blanket form takes theirs down with yours.
`ps -eo pid,etime,command | grep magma.exe` first, then `kill <pid>`.

**`import` defeats `AttachSpec`'s laziness.** `AttachSpec` loads packages on demand, but
`import "X.m" : f;` compiles `X.m` immediately as its own package — so intrinsics from *other*
spec files are unresolved inside it and you get `Undefined reference` at call time. Touch one
intrinsic from the needed package first:

    AttachSpec("ShimuraQuotients.spec");
    _ := ClassNumberLU(-4);              // forces ClassNumberData.m to load
    import "TraceFormula.m" : Hurwitz;   // now resolves

## Normaliz (required for any polytope solve)

The polymake backend is dead; `nmzsolve.py` replaced it. `NORMALIZ_BIN` **must** be set — the
fallback path it computes does not exist in this checkout:

    export NORMALIZ_BIN=~/Documents/GitHub/normaliz-3.11.1/normaliz

Without it, or above the cached frontier, a fresh solve fails *silently* — you get "no solutions"
rather than an error, and a partially-cached base returns a wrong answer instead of complaining.
Committed cache files use escaped `\[` line starts, so byte-comparing them against fresh Normaliz
output fails on identical content; compare as vector sets.

## Where things live

    main                the working branch: code + paper/ + PLAN.md/HANDOFF.md
    m0-theta-campaign   research data, probes and triage tooling (vvdata/weyl-campaign/)

Worktrees live under `worktrees/` inside this repo checkout, not as siblings of it — deliberate,
to keep the `GitHub/` directory tidy:

    .                    main (this checkout)
    worktrees/campaign   m0-theta-campaign

**When adding a new worktree, put it under `worktrees/<name>` here**, e.g.:

    git worktree add worktrees/<name> <branch>

**Two long-lived branches: `main` and `m0-theta-campaign`**, plus short-lived topic branches
while their PRs are open. Everything else is retired and
preserved as an `archive/<name>` tag on `origin` (10 of them). `whbasis-speedup` went too — its
one commit was already in `main` via the `04f1d7b` cherry-pick.

**⚠ The campaign branch carries a FULL CODE TREE, not just data.** So a probe run from
`worktrees/campaign` uses *that branch's* code, not `main`'s. It had drifted 103 commits behind
before being merged up on 2026-09-04; **merge `main` into it before trusting any measurement
taken there.**

**The invariant that keeps this from biting again — check it, don't rely on discipline.** The
gaps above happened in files at **SHARED PATHS**: paths that exist on *both* branches and can
therefore drift apart silently (`nmzsolve.py` at the root, `vvdata/gtsweep.m`). Anything under
`vvdata/weyl-campaign/` can never diverge, because `main` does not have it. So:

    git diff origin/main origin/m0-theta-campaign --name-only -- ':!vvdata/weyl-campaign/*'

**should print nothing but doc files.** Anything else is a silent divergence — run it before
trusting either branch's code. Had this existed, the nine-day `nmzsolve.py` gap would have shown
up immediately. Corollary: **make a change to a shared-path file via a PR to `main`, then merge `main` down.**
If something belongs only to the research line, put it under `vvdata/weyl-campaign/` — that is
why the FIRE variant is `vvdata/weyl-campaign/gtsweep_fire.m` and not a fork of
`vvdata/gtsweep.m`.

`nmzsolve.py` used to conflict on that merge, because campaign carried the **t-shift fallback**
and `main` did not. **Resolved 2026-09-04: the fallback is now IN `main`** (with the two files it
reads at runtime, `polymake/tshift_{core,w0}_420.txt`); the two copies are identical, so that
conflict should not recur. Its generator and probe stay on campaign
(`vvdata/weyl-campaign/tshift_gen.py`, `tshift308.m`) per the tooling convention. The fallback is
guarded — it needs `m_pole==0 && k24==12 && cuspidal==0`, a `polymake/tshift_w0_<M>.txt`, and a
lower cached rung — so it is inert at levels without those files, and it is validated at
**M = 420 only**.

**⚠ `tier1-models` is RETIRED (2026-09-04) and `main` carries everything it had.** It was
fast-forwarded into `main` — the two were the same commit — and then deleted, local and remote,
along with the now-redundant `worktrees/mainport`. Older material (`HANDOFF.md`, `PLAN.md`,
memory) still says things like "`main` is code only", "the `tier1-models` merge trap", or refers
to `-campaign` / `-mainport` sibling directories. **All of that is historical** — there is one
code branch now, and it is `main`. Do not recreate `tier1-models`; work for `main` goes on a
topic branch and lands by PR.

**Triage tooling lives on the campaign branch under `vvdata/weyl-campaign/`, never at the repo
root.** So `git log --all -- cmsupply.m` reports "not in any branch" for a file that is committed.
Search by basename:

    git log --all --oneline --name-only --pretty=format: -- '*cmsupply.m'

Before claiming a file exists nowhere in git, search that way — a root-path query is guaranteed to
miss it.

## Conventions

* Scratch scripts belong in `vvdata/weyl-campaign/` on the campaign branch, not `/tmp` — `/tmp` is
  purged nightly and has already eaten one driver that had to be rewritten from a handoff.
* Record *why* a probe's number is trustworthy, next to the probe. Several results here are proxies
  with a limited domain of validity, and the failure mode is quoting one outside it.
