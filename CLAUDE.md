# Working in this repo

`README.md` covers loading the curve data. This file covers the things that are true regardless of
what you are working on, and that have each cost a session at least once.

## Start here

* **`PLAN.md`** — what to do next. Opens with a "Picking this up cold" block.
* **`HANDOFF.md`** — what happened. Authoritative on state; it wins over `PLAN.md` and over any
  agent memory when they disagree.

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

**Two branches, that is all: `main` and `m0-theta-campaign`.** Everything else is retired and
preserved as an `archive/<name>` tag on `origin` (10 of them). `whbasis-speedup` went too — its
one commit was already in `main` via the `04f1d7b` cherry-pick.

**⚠ The campaign branch carries a FULL CODE TREE, not just data.** So a probe run from
`worktrees/campaign` uses *that branch's* code, not `main`'s. It had drifted 103 commits behind
before being merged up on 2026-09-04; **merge `main` into it before trusting any measurement
taken there.** One conflict recurs: `nmzsolve.py`, where the campaign copy is a strict superset
(the t-shift fallback plus its `tshift` helper files, which `main` has never received) — keep the
campaign side.

**⚠ `tier1-models` is RETIRED (2026-09-04) and `main` carries everything it had.** It was
fast-forwarded into `main` — the two were the same commit — and then deleted, local and remote,
along with the now-redundant `worktrees/mainport`. Older material (`HANDOFF.md`, `PLAN.md`,
memory) still says things like "`main` is code only", "the `tier1-models` merge trap", or refers
to `-campaign` / `-mainport` sibling directories. **All of that is historical** — there is one
code branch now, and it is `main`. Do not recreate `tier1-models`; commit to `main` directly.

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
