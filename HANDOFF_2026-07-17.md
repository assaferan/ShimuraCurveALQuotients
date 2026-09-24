# Handoff — producing more X_0(D,N)* cover models

**Date:** 2026-07-17. **Written for:** the next agent, on a *different machine* (so nothing in
`/tmp` or the scratchpad transfers — everything you need is in this file or in git).

**Overarching goal:** produce more subhyperelliptic cover models for Shimura-curve Atkin-Lehner
star quotients X_0(D,N)*, written to `data/models/models_D_N.m`.

> (This supersedes an older HANDOFF.md about the Weil-budget / filtering-pipeline task, which is
> long finished. If you need that history it's in git log on `main`/earlier branches.)

---

## 0. IMMEDIATE STATE / ACTION NEEDED

- **Two branches matter:**
  - `special-fiber-modular` — the mainline. All this session's *finished* work is committed AND
    pushed here (last commit `4e6f93b "adding model verification to CI"`). There is an open PR to
    `main`; CI runs on every push.
  - `more-models` — a **LOCAL, UNPUSHED** branch off `special-fiber-modular` where the "produce more
    models" investigation lives. Nothing is committed on it yet (work so far = analysis + two
    untracked model files). **Push it from the old machine if you want the branch, or just recreate
    it from `special-fiber-modular` — no commits are lost either way.**
- **Two new models are UNTRACKED and only on the old machine's disk:**
  `data/models/models_6_5.m` (15 cover-keys) and `data/models/models_14_5.m` (14 cover-keys),
  produced with NO code change. If they weren't committed/pushed, **regenerate them** (§4) — cheap.
- Also untracked (pre-existing, leave alone): `data/cn16_discriminants.txt`, one modified
  `polymake/polymake_script_420_145_0` (a benign cache write).

---

## 1. WHAT WAS ACCOMPLISHED THIS SESSION (committed on special-fiber-modular)

All DONE, validated, committed, pushed:

1. **CM-fetch speedup** (`6e2eab8`) + **Prop 5.6 existence fix** (`f32c0c7`): ~210x speedup of
   `FieldsOfDefinitionOfCMPointFast` (field-of-def degree read from group orders, no number field
   built) plus a `GCD(D,f)=1` phantom-CM-point guard.

2. **Solve-stage sign disambiguation** (`a919182`): the "not-enough-CM" blocker was mis-diagnosed.
   Real cause: `find_hauptmodul_signs_quadratic` (SchoferFormula.m) filtered the star hauptmodul
   s(P)'s minpoly against a COVER's field of definition, which is often **degree 4** — far too loose,
   so spurious sign candidates survived, the point was discarded, and the quadratic supply was
   "exhausted" ("No possible choices of CM points left"). Fix:
   - Filter against the **star curve's** degree-2 field of definition (where s(P) actually lives).
   - Genuine degree-2 ambiguities (both sign choices have a root in K) are recorded on the table
     (`AmbiguousSigns` attr, cmtables.m) and resolved at the **solve stage**:
     `AllEquationsAboveCovers` searches the small product of per-point sign choices and keeps the
     assignment determining the most covers **of correct genus** (a wrong sign makes a quadratic
     constraint inconsistent, deferring that cover).

3. **Direct low-genus quotient** (`a919182`): `backfill_deferred` used
   `CurveQuotient(AutomorphismGroup(C,[w]))`, which Magma rejects for genus <= 1 over Q. Added
   `direct_involution_quotient(C, w)` (EquationsCovers.m): for C: y^2=f(x) and coordinate involution
   w:(x,y)->(mu(x),rho(x)y), build the quotient from invariants (U = non-constant symmetric function
   of {x,mu(x)}, V=(1+rho)y, V^2=G(U)). Verified isomorphic to the AutomorphismGroup result on
   y^2=x^6+1. Kept the AutGrp path for genus > 1.

4. **Model verification** (`784c8cb`, `4e6f93b`): `ModelVerification.m` + `tests/ModelChecks.m` (CI
   target). Four checks, none using the Borcherds/Schofer machinery that built the models:
   [1] genus self-consistency, [2] genus vs Shimura genus formula, [3] Weil-polynomial divisibility
   across nested cover keys, [4] **trace-formula point count**:
   `#(X/W)(F_p) == ComputePointsViaTrace(X,p,1)` (Eichler-Selberg via `TraceDNewALFixed`, NO modular
   symbols). Result: 71 files, 7087 checks, 0 failures, ~6 min; a single perturbed coefficient is
   rejected (13 failures).
   - **Magma quirks learned (documented in the test):** `load` needs a literal filename (model files
     read via `Read`+`eval`); an `eval` inside a procedure that closes over an outer variable
     **segfaults Magma 2.29** (keep test top-level); a `{...}` docstring cannot contain nested
     braces; `ModularSymbols` spaces are cached forever and exhaust memory across many levels (hence
     check [4] uses the trace formula).

**5 bases that previously produced NOTHING now produce verified models:**
22_7 (15 keys), 22_13 (11), 38_7 (8), 14_3 (6), 34_5 (6).
**Still timing out** (highest demand): **14_37, 34_11** — adequate CM supply, just need a long
uninterrupted run on a free machine.

---

## 2. THE STRATEGIC FINDING (the meat of the next step)

A coverage audit (script in §5) found:

- **423** genuine targets = D>1 star curves with >=1 subhyperelliptic cover.
- **71** have a model file (17%). **352 are MISSING**, by demand `num_vals = max(2*g+5)`:

  | demand | missing | | demand | missing |
  |--------|---------|-|--------|---------|
  | 7      | 16      | | 13     | 47      |
  | **9**  | **173** | | 15     | 34      |
  | 11     | 54      | | 17+    | 28      |

  **189 of the 352 sit at demand 7-9** — the regime the current machinery handles.

- **Existing model sets are typically PARTIAL** (empty cover-keys = deferred, unrecovered), e.g.
  committed `10_41`: 7/14 empty, `14_29`: 7/14. Sparse output is in-family, not a defect. "More
  models" = both more bases AND filling in existing sets.

### Decisive experiment: the missing bases FAIL (not never-run), in DIVERSE ways.

Ran the pipeline on ~23 of the cheapest low-demand missing bases. Of the 13 that finished cleanly:

- **2 SUCCEEDED with no code change**: `6_5` (15 keys), `14_5` (14 keys).
- The rest failed with **>=5 distinct errors** — a long tail of bugs, NOT one wall:

  | failure | count | example bases |
  |---------|-------|---------------|
  | `Runtime error in assert: Assertion failed` | 3 | 10_9, 15_4, 21_4 |
  | `Could not find enough rational CM points!` (in BorcherdsForms, UPSTREAM) | 2 | 22_3, 34_3 |
  | `RationalNumber: s does not represent a rational number!` | 2 | 21_1, 51_1 |
  | `Could not find enough points, sorry!` | 1 | 26_3 |
  | `^: Argument 2 is too large` | 1 | 33_1 |

  **Histogram caveat:** 12 of the 23 probes were killed by a 3600s timeout with **0-byte output**
  (genmodels writes only at the very end, so a killed run logs nothing). Those are UNMEASURED, not
  failures. A verbose rerun (`SetVerbose("ShimuraQuotients",1)` so errors stream even if killed) was
  in flight when this machine was vacated — **redo it** to get the true ranking before choosing a fix.

---

## 3. RECOMMENDED NEXT STEP

1. **Rebuild the failure histogram cleanly** on a free machine: run the ~30 cheapest missing bases
   (§6) with verbose logging + a generous timeout; count error classes. (The 12 unmeasured could
   shift the ranking.)
2. **Fix the single largest fixable class.** Front-runner: bare **`Assertion failed`** (a plain
   assert localizes to one line — cheapest to diagnose). Get the line via `SetDebugOnError(true)` on
   `10_9` (see the debug driver in §4).
3. Re-run affected bases; measure unlocked count; iterate to the next class.
4. `Could not find enough rational CM points` is upstream in `BorcherdsForms` and may be its own
   project (see memory `cm-point-supply-lever`, `fieldsofdef-maxdegree-fetch-speedup`).
5. When a class works, **run at scale** over the 189 demand-7/9 targets and **verify every output**
   with `tests/ModelChecks.m` (it auto-discovers new model files) before trusting them.

**Do NOT** assume this session's "No possible choices" fix applies — that's DONE; the remaining
failures are different, earlier errors.

---

## 4. HOW TO RUN THINGS (environment)

- **Magma via sage:** `sage -sh -c "magma -b <args> file.m"` (also `/opt/magma/current/magma`).
  Batch `-b`; params as `name:=value` read with `StringToInteger(...)`.
- **Compile check:** `AttachSpec("ShimuraQuotients.spec"); printf "SPEC-OK\n"; quit;`
- **Run a CI test target:** `magma -b verbose:=3 target:=ModelChecks.m exitsignal:="" run_tests.m`
- **Shared machine reality:** 256 cores but often saturated by another user's ~200 jobs. Run at
  **normal priority** (`nice -18` starved us badly; single-core jobs are negligible). Don't trust
  wall-clock under load. genmodels writes model output only at the END → a killed run leaves an empty
  file; use verbose logging to capture progress.

### The genmodels harness (RECREATE — it lived at `/tmp/genmodels.m`, NOT in git)

Run from repo root: `sage -sh -c "magma -b D_s:=<D> N_s:=<N> /tmp/genmodels.m"`

```magma
SetQuitOnError(true);
SetColumns(0);
AttachSpec("ShimuraQuotients.spec");
Px<x> := PolynomialRing(Rationals());
D := StringToInteger(D_s); N := StringToInteger(N_s);
// Skip known-deferred vx bases (memory vx-laurent-n0-circular): huge Zero-side space then
// crash on the AbsEltseq "vx ge 0" assert. Fail fast.
vx_skip := {<95,1>, <115,1>, <123,1>, <129,1>};
if <D,N> in vx_skip then WriteStderr(Sprintf("SKIP vx base %o_%o\n", D, N)); quit; end if;
curves := GetHyperellipticCandidates();
Xstar := rep{X : X in curves | X`D eq D and X`N eq N and IsStarCurve(X)};
t0 := Realtime();
covers, ws := AllEquationsAboveCovers(Xstar, curves : Prec := 100);
WriteStderr(Sprintf("computed AllEquationsAboveCovers in %os\n", Realtime()-t0));
agg := AssociativeArray();
for label in Keys(covers) do
  X := curves[label]; Wkey := Sort([w : w in X`W]);
  if not IsDefined(agg, Wkey) then agg[Wkey] := [* *]; end if;
  for base in Keys(covers[label]) do Append(~agg[Wkey], <X`g, covers[label][base]>); end for;
end for;
keystr := func<W | "[Integers()|" cat (#W eq 0 select "" else &cat[IntegerToString(w) cat (i lt #W select "," else "") : i->w in W]) cat "]">;
header := "models := AssociativeArray();";
lines := [
  Sprintf("// Subhyperelliptic cover models for X_0(%o,%o)*", D, N),
  "// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).",
  "P<x> := PolynomialRing(Rationals());",
  header
];
ncov := 0;
for Wkey in Keys(agg) do
  entries := [];
  for m in agg[Wkey] do
    g := m[1]; C := m[2];
    if Type(C) eq CrvHyp then
      f, h := HyperellipticPolynomials(C);
      Append(~entries, Sprintf("<%o, P!%o, P!%o>", g, Coefficients(f), Coefficients(h)));
    else
      Append(~entries, Sprintf("<%o, \"CRV\", %m>", g, [Sprint(p) : p in DefiningPolynomials(C)]));
    end if;
  end for;
  Append(~lines, Sprintf("models[%o] := [* %o *];", keystr(Wkey), Join(entries, ", ")));
  ncov +:= 1;
end for;
fname := Sprintf("data/models/models_%o_%o.m", D, N);
Write(fname, Join(lines, "\n") : Overwrite);
WriteStderr(Sprintf("WROTE %o (%o cover-keys)\n", fname, ncov));
quit;
```

### Debug driver for a FAILING base (get error + stack)

```magma
SetColumns(0); AttachSpec("ShimuraQuotients.spec");
SetVerbose("ShimuraQuotients",1); SetDebugOnError(true);
curves := GetHyperellipticCandidates();
X := rep{Y : Y in curves | Y`D eq 10 and Y`N eq 9 and IsStarCurve(Y)};
covers, ws := AllEquationsAboveCovers(X, curves : Prec := 100);
```

---

## 5. THE AUDIT SCRIPT (RECREATE — lived in scratchpad, NOT in git)

Lists missing subhyperelliptic targets by demand. Run from repo root; writes `MISSING_TARGETS.txt`
with lines `base demand level`:

```magma
SetColumns(0); AttachSpec("ShimuraQuotients.spec");
curves := GetHyperellipticCandidates();
stars := [X : X in curves | IsStarCurve(X) and X`D gt 1];
have := Split(Pipe("ls data/models/models_*.m 2>/dev/null", ""), "\n");
have := [f : f in have | #f gt 0];
havekeys := {};
for f in have do
  parts := Split(f, "/"); nm := parts[#parts]; core := nm[8..#nm-2]; dn := Split(core, "_");
  if #dn eq 2 then Include(~havekeys, <StringToInteger(dn[1]), StringToInteger(dn[2])>); end if;
end for;
out := "";
for X in stars do
  subcovers := [i : i in X`CoveredBy | assigned curves[i]`IsSubhyp and curves[i]`IsSubhyp];
  if #subcovers eq 0 then continue; end if;
  if <X`D, X`N> in havekeys then continue; end if;
  gl := [curves[i]`g : i in subcovers]; nv := Maximum([2*g+5 : g in gl]);
  out cat:= Sprintf("%o_%o %o %o\n", X`D, X`N, nv, X`D*X`N);
end for;
Write("MISSING_TARGETS.txt", out : Overwrite); quit;
```
The `IsSubhyp` attribute (set by `UpdateByGenus` in ShimuraQuotients.m) is what turns "1945 D>1
stars" into "423 genuine targets".

---

## 6. TARGET LISTS

**All 16 demand-7 missing bases** (cheapest; * = already produces a model):
```
6_5* 6_25 6_35 10_9 10_21 14_5* 15_2 21_1 21_2 22_3 33_1 33_2 34_3 57_5 210_1 330_1
```

**Cheapest 30 demand<=9 by level D*N** (good first batch to run + histogram):
```
21_1 15_2 6_5 33_1 21_2 51_1 57_1 15_4 22_3 33_2 69_1 14_5 35_2 26_3 39_2 21_4 85_1
10_9 91_1 34_3 51_2 15_7 21_5 22_5 38_3 57_2 123_1 14_9 65_2
```

**Known error mapping from the probe** (start debugging here):
- `Assertion failed`: 10_9, 15_4, 21_4
- `Could not find enough rational CM points` (BorcherdsForms, upstream): 22_3, 34_3
- `s does not represent a rational number`: 21_1, 51_1
- `Could not find enough points, sorry!`: 26_3
- `^: Argument 2 is too large`: 33_1

---

## 7. CONVENTIONS & MEMORY

- **Git:** commit messages have NO `Co-Authored-By` trailer. Assistant typically STAGES; USER
  commits/pushes (controls CI pushes). Models committed locally by the user.
- **Polymake cache:** `polymake/` has ~732 tracked cache files. Untracked ones split into "has a
  tracked script<->solution mate" (legit, keep/stage) vs "orphans" (BorcherdsForms runaway-bug
  artifacts, up to MBs, safe to `git clean`). Memory `polymake-cache-keep`.
- **Persistent memory:**
  `/home/assaferan/.claude/projects/-home-assaferan-GitHub-ShimuraCurveALQuotients/memory/`
  (MEMORY.md index). Relevant: `fieldsofdef-maxdegree-fetch-speedup`, `cm-point-supply-lever`,
  `vx-laurent-n0-circular`, `polymake-cache-keep`, `special-fiber-modular-branch`. Verify
  file/function names still exist before relying on them.
- **Files changed this session** (committed on special-fiber-modular): `cmtables.m`,
  `SchoferFormula.m`, `EquationsCovers.m`, `ModelVerification.m`, `ShimuraQuotients.spec`,
  `tests/ModelChecks.m`.
- **`pkill` self-matches its own command line** (exit 143/144) — use exact patterns or pid loops.

---

## 8. ONE-LINE SUMMARY

Session delivered sign-disambiguation + low-genus-quotient fixes (5 new verified models) and a
4-check model verifier in CI. Next: the audit shows **352 missing targets, 189 at demand 7-9**,
failing in a **long tail of ~5+ distinct bugs** (2 bases already work untouched). Rebuild the
failure histogram cleanly, fix the largest class (likely the bare `Assertion failed`), and scale.
