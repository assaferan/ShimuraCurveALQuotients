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

* **`HANDOFF.md`** — what happened. Authoritative on state; it wins over any agent memory when
  they disagree.

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
* Delete a topic branch once its PR merges.

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
2026-09-23: both runs were killed by macOS for memory pressure at exactly `X0_206_1.m`. So a local
`run_tests.m` with no target learns nothing past the `X0_1*` range. Either use `target:=` /
`filename:=` for the tests you care about, or run it on a server. ⚠ A killed suite LOOKS like a
clean one — check that every test file printed a result line. `run_tests.m` prints `Tests failed:`
only when something failed, so its absence alone proves nothing.

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

The polymake backend is dead; `nmzsolve.py` replaced it. `NORMALIZ_BIN` **must** point to a
Normaliz binary — the fallback path it computes does not exist in this checkout. On Eran's machine:

    export NORMALIZ_BIN=~/Documents/GitHub/normaliz-3.11.1/normaliz

Elsewhere, find the local install. CI uses `/usr/bin/normaliz` from the apt package `normaliz-bin`.

Without it, or above the cached frontier, a fresh solve fails *silently* — you get "no solutions"
rather than an error, and a partially-cached base returns a wrong answer instead of complaining.
Committed cache files use escaped `\[` line starts, so byte-comparing them against fresh Normaliz
output fails on identical content; compare as vector sets.

## Where things live

**Don't assume which branches exist; check `git branch -r`.** `main` is the default PR target,
but some PRs target other branches (for example `integration`, which stages fixes for a combined
pipeline rerun), so read a PR's base before reasoning about it. Work goes on a topic branch and
lands by PR. Retired branches are preserved as `archive/<name>` tags on `origin`.

`git worktree list` shows what is checked out where. New worktrees go under `worktrees/<name>`
in this checkout. For the campaign branch (research data, probes and triage tooling under
`vvdata/weyl-campaign/`):

    git worktree add worktrees/campaign m0-theta-campaign

**Never `git pull` a clone with jobs running from it.** `AttachSpec` loads packages as they are
first used, so updating the tree under a running job mixes two versions of the code.

**⚠ The campaign branch carries a FULL CODE TREE, not just data.** So a probe run from
`worktrees/campaign` uses *that branch's* code, not `main`'s. **Merge `main` into it before
trusting any measurement taken there.**

**The invariant that keeps this from biting again — check it, don't rely on discipline.** The
gaps above happened in files at **SHARED PATHS**: paths that exist on *both* branches and can
therefore drift apart silently (`nmzsolve.py` at the root, `vvdata/gtsweep.m`). Anything under
`vvdata/weyl-campaign/` can never diverge, because `main` does not have it. So:

    git diff origin/main origin/m0-theta-campaign --name-only -- ':!vvdata/weyl-campaign/*'

**should print nothing but doc files.** Anything else is a silent divergence — run it before
trusting either branch's code. Corollary: **make a change to a shared-path file via a PR to `main`, then merge `main` down.**
If something belongs only to the research line, put it under `vvdata/weyl-campaign/` — that is
why the FIRE variant is `vvdata/weyl-campaign/gtsweep_fire.m` and not a fork of
`vvdata/gtsweep.m`.

`nmzsolve.py` on `main` carries the **t-shift fallback**, with the two files it reads at runtime,
`polymake/tshift_{core,w0}_420.txt`. Its generator and probe stay on campaign
(`vvdata/weyl-campaign/tshift_gen.py`, `tshift308.m`) per the tooling convention. The fallback is
guarded — it needs `m_pole==0 && k24==12 && cuspidal==0`, a `polymake/tshift_w0_<M>.txt`, and a
lower cached rung — so it is inert at levels without those files, and it is validated at
**M = 420 only**.

`tier1-models` is retired and merged into `main`; do not recreate it. Older notes mentioning it or
`-campaign` / `-mainport` sibling directories are historical.

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
* **Every number in a test states its provenance**: the paper with table or equation, an LMFDB
  label, or the independent command that produced it (another Magma method, Sage/PARI, someone
  else's code). "Computed by this repo at commit X" makes a regression pin, not a verification, so
  label it that way.
* **Don't edit a failing test to make it pass.** A failing test is evidence. Fix the code, or, if
  you think the test itself is wrong, stop and raise it with a person, giving the independent
  source that shows the expected value is wrong. Never change an expected value, loosen an assert,
  or delete a check on your own.

## Writing PRs, comments and replies

* **PR bodies above all: write for a collaborator reading cold, not for another agent.** Open
  with two or three plain sentences saying what was wrong and what the change does, the
  mathematical point before the implementation detail. Then one bullet per file and a line or two
  on how it was verified. Length follows importance: a small fix needs a few lines, while an
  important problem can take as long as it needs, provided a reader can follow it. Leave out
  internal names, offsets, commit-by-commit history and evidence the reader doesn't need.
* The same plain style applies elsewhere, more briefly:
  * Code comments: 1–4 lines stating the constraint the code cannot show.
  * Tracker comments: a few plain lines on what changed, not a report.
* **Don't comment idiomatic patterns.** A comment must state a constraint that the code cannot
  show and that a competent reader of this codebase would not already know. When unsure, leave it
  out and let the reviewer ask.
* **Describe an implementation on its own terms**, not by criticising an alternative. Say what it
  does, how to call it, and what it costs. Mention comparable APIs as neighbours, not as warnings.
