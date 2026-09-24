# Target-restricted ratpts sweeps — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the Table 1/6/10 rational-point sweeps cheaper and more robust by (a) restricting `EquationsOfCovers` to the target covers (lowers the CM-point demand), (b) a cheap CM-point pre-check that skips genuinely-unsolvable groups fast, and (c) a portable per-group OS timeout so a single hang can't abort a sweep.

**Architecture:** Three independent layers over the existing Magma pipeline. The math change is a single new optional `Targets` parameter on `EquationsOfCovers(Xstar, curves)` that (1) restricts `num_vals`/`MaxNum` to the target immediate covers and (2) restricts the per-cover solve by filtering `SchoferTable`K_idxs`. A new `EnoughCMPointsForTargets` intrinsic gives the cheap guard. A `run_table.sh` shell runner wraps each group in a `perl alarm` timeout (no `timeout` binary on this box) with resume + TIMEOUT logging. Drivers (`ratpts_table1/6/10.m`) pass row targets, call the guard, and log distinct verdicts.

**Tech Stack:** Magma (`.m`, AttachSpec `ShimuraQuotients.spec`); CHIMP (Table 10 only); bash + perl for the runner. No pytest — verification is via small `magma` invocations whose printed output is the oracle. The settled cases `D=6,N=29` (BIELLIPTIC) and `D=34,N=5` (insufficient CM, 4<9) are the test fixtures.

**Key code facts (from reading the codebase + `spike_lattice.m`):**
- `curves` is indexed by `CurveID`; `Xstar`CoveredBy` and `SchoferTable`Keys_fs` (positive entries) are `CurveID`s that index `curves`.
- `EquationsCovers.m:178` sets `num_vals := Maximum([2*g+5 : g in [curves[i]`g : i in Xstar`CoveredBy]])` — the global demand.
- `cmtables.m:22` sets `K_idxs := [i : i->k in keys_fs | k gt 0]`; the solve loops (`RationalConstraintsOnEquations` EquationsCovers.m:8, the kernels loop EquationsCovers.m:131) iterate `K_idxs`. `RationalConstraintsOnEquations` has `require #ds ge 2*g+3` (EquationsCovers.m:14) per cover, so non-target high-genus siblings MUST be removed from `K_idxs` before that call.
- `RationalandQuadraticCMPoints(X : bd:=8, coprime_to_level:=true)` (ShimuraQuotients.m:1403) returns `rat, quad`; ~27 s for `D=34,N=5`, requires `X` to be a star quotient.
- No `timeout`/`gtimeout` installed; `perl` is available.

---

## Phase A — Timeout & runner infrastructure (independent, low-risk)

### Task A1: Ordered `TABLE6` list + `idx:=` support in `ratpts_table6.m`

**Files:**
- Modify: `ratpts_table6.m` (currently selects via `CANDIDATES`; `ratpts_table1/10.m` already accept `idx:=`)

- [ ] **Step 1: Add an ordered `TABLE6` and an `idx:=` entry path.** After the existing `CANDIDATES` block, add a flat ordered list mirroring the `CANDIDATES` rows (one `<D,N,gensets>` per group, D*N ascending) and an entry point so a single group can be run by index:

```
// Ordered, indexable view of the Table 6 groups (D*N ascending), for run_table.sh.
TABLE6 := CANDIDATES;  // CANDIDATES is already the ordered <D,N,[gensets]> list
Sort(~TABLE6, func< a, b | (a[1]*a[2]) ne (b[1]*b[2]) select (a[1]*a[2])-(b[1]*b[2])
                           else (a[2] ne b[2] select a[2]-b[2] else a[1]-b[1]) >);

if assigned idx then
    i := StringToInteger(idx);
    printf "Running single TABLE6 group #%o of %o.\n", i, #TABLE6;
    CANDIDATES := [* TABLE6[i] *];   // restrict the main loop to this one group
end if;
```

(The existing `for entry in CANDIDATES` loop then runs just that group when `idx` is set, or all groups otherwise — no other change to the loop body.)

- [ ] **Step 2: Verify single-group selection works (no Borcherds compute needed to check indexing).**

Run: `magma idx:=1 ratpts_table6.m 2>&1 | grep -E "Running single TABLE6|==== D="`
Expected: prints `Running single TABLE6 group #1 of N.` followed by exactly one `==== D=.. N=.. ====` header.

- [ ] **Step 3: Commit.**

```bash
git add ratpts_table6.m
git commit -m "ratpts_table6: add ordered TABLE6 + idx:= single-group selection"
```

### Task A2: `run_table.sh` — portable per-group timeout, resume, TIMEOUT logging

**Files:**
- Create: `run_table.sh`

- [ ] **Step 1: Write the runner.** It loops group indices, wraps each `magma idx:=i` in a `perl alarm` timeout, tee's output, and records TIMEOUT rows. Resume is by checking the results file for an already-present `(D,N)` for that index (parsed from a dry `magma idx:=i` header is overkill; instead the driver itself already skips settled `(D,N,W)` rows, so the runner only needs to record TIMEOUTs and continue).

```bash
#!/usr/bin/env bash
# Usage: run_table.sh <table:1|6|10> <timeout-seconds> [lo] [hi]
# Runs each group as its own `magma idx:=i ratpts_table<table>.m`, killed after
# <timeout-seconds> via perl alarm (no coreutils `timeout` on this box). A hang/OOM
# kills only that process; the sweep continues. Per-group stdout -> table<table>_run_*.log
set -u
TBL="$1"; T="$2"; LO="${3:-1}"; HI="${4:-}"
SCRIPT="ratpts_table${TBL}.m"
RESULTS="ratpts_table${TBL}_results.txt"
STAMP="$(date +%Y%m%d_%H%M%S)"
LOG="table${TBL}_run_${STAMP}.log"
[ -z "$HI" ] && HI=$(magma -b show:=1 "$SCRIPT" 2>/dev/null | grep -Eo 'groups' >/dev/null; echo "")
# If HI not given, the driver prints its group count; fall back to a large cap and let
# out-of-range idx error per group (cheap) — simpler: require HI for table 6/1, optional for 10.
if [ -z "$HI" ]; then echo "Provide hi (group count) explicitly."; exit 2; fi
echo "Sweeping table $TBL groups $LO..$HI, timeout ${T}s each. Log: $LOG"
for i in $(seq "$LO" "$HI"); do
  echo "==== group #$i (timeout ${T}s) ====" | tee -a "$LOG"
  perl -e 'alarm shift; exec @ARGV' "$T" magma idx:="$i" "$SCRIPT" >>"$LOG" 2>&1
  rc=$?
  if [ "$rc" -eq 142 ] || [ "$rc" -eq 124 ]; then   # 142 = 128+SIGALRM(14)
    echo "  group #$i TIMEOUT after ${T}s (rc=$rc)" | tee -a "$LOG"
    printf '# TIMEOUT\tgroup=%s\ttable=%s\t>%ss\n' "$i" "$TBL" "$T" >> "$RESULTS"
  fi
done
echo "Sweep done. See $LOG and $RESULTS"
```

- [ ] **Step 2: Make executable and smoke-test the timeout path with a trivial sleep.**

```bash
chmod +x run_table.sh
perl -e 'alarm shift; exec @ARGV' 2 sleep 30; echo "rc=$?"
```
Expected: returns after ~2 s with a non-zero `rc` (the `alarm` killed `sleep`), confirming the wrapper kills a long child.

- [ ] **Step 3: Commit.**

```bash
git add run_table.sh
git commit -m "Add run_table.sh: portable perl-alarm per-group timeout sweep runner"
```

---

## Phase B — Cheap CM-point pre-check

### Task B1: `EnoughCMPointsForTargets` intrinsic

**Files:**
- Modify: `ShimuraQuotients.m` (add intrinsic immediately after `RationalandQuadraticCMPoints`, ~line 1462)

- [ ] **Step 1: Add the intrinsic.**

```
intrinsic EnoughCMPointsForTargets(Xstar::ShimuraQuot, curves::SeqEnum[ShimuraQuot], Targets::SetEnum)
          -> BoolElt, RngIntElt, RngIntElt
{Cheap predicate: are there enough rational+quadratic CM points to solve the target immediate
 covers? Returns (available ge required), required, available. required = max(2g+5) over the
 immediate covers of Xstar whose W is in Targets.}
    target_keys := [i : i in Xstar`CoveredBy | curves[i]`W in Targets];
    require #target_keys gt 0 : "No immediate cover of Xstar matches Targets.";
    required := Maximum([2*curves[i]`g + 5 : i in target_keys]);
    rat, quad := RationalandQuadraticCMPoints(Xstar : bd := 8, coprime_to_level := true);
    available := #rat + #quad;
    return available ge required, required, available;
end intrinsic;
```

- [ ] **Step 2: Verify on the two fixtures.** Write `verify_guard.m`:

```
AttachSpec("ShimuraQuotients.spec");
SetVerbose("ShimuraQuotients", 0);
curves := GetHyperellipticCandidates();
procedure chk(D, N, gens, expect_ok)
    Xstar := rep{X : X in curves | X`D eq D and X`N eq N and IsStarCurve(X)};
    T := { AllALsFromGens(gens, D*N) };
    ok, req, avail := EnoughCMPointsForTargets(Xstar, curves, T);
    printf "D=%o N=%o gens=%o : ok=%o (req=%o avail=%o) expected_ok=%o  %o\n",
        D, N, gens, ok, req, avail, expect_ok, (ok eq expect_ok select "PASS" else "FAIL");
end procedure;
chk(34, 5, {2,5}, false);    // 4 < 9  -> ok=false
chk(10, 7, {2,5}, true);     // genus-0 target, plenty -> ok=true
exit;
```

Run: `magma verify_guard.m 2>&1 | grep -E "PASS|FAIL"`
Expected: both lines `PASS`; the `34,5` line shows `req=9 avail=4`.

- [ ] **Step 3: Commit.**

```bash
git add ShimuraQuotients.m verify_guard.m
git commit -m "Add EnoughCMPointsForTargets cheap CM-point pre-check + verify script"
```

---

## Phase C — Target-restricted construction in `EquationsOfCovers`

### Task C1: `Targets` parameter restricting `num_vals` and the solve set

**Files:**
- Modify: `EquationsCovers.m:173-190` (the `EquationsOfCovers(Xstar, curves : Prec)` intrinsic)

- [ ] **Step 1: Add the `Targets` parameter and restrict `num_vals` + `K_idxs`.** Replace the body so that, when `Targets` is non-empty, the genus list comes from the target covers and the solved `K_idxs` is filtered to target covers before the solve:

```
intrinsic EquationsOfCovers(Xstar::ShimuraQuot, curves::SeqEnum[ShimuraQuot]
                            : Prec := 100, Targets := {}) -> SeqEnum, Assoc, SeqEnum
{Determine the equations of the immediate covers of X. If Targets (a set of W subgroups) is
 nonempty, only those covers are solved and the CM-point demand is set from their genera.}
    fs := BorcherdsForms(Xstar, curves : Prec := Prec);
    d_divs := &cat[[T[1]: T in DivisorOfBorcherdsForm(f, Xstar)] : f in [fs[-1], fs[-2]]];
    all_cm_pts := CandidateDiscriminants(Xstar, curves);

    if #Targets gt 0 then
        cover_keys := [i : i in Xstar`CoveredBy | curves[i]`W in Targets];
        require #cover_keys gt 0 : "Targets match no immediate cover of Xstar.";
    else
        cover_keys := [i : i in Xstar`CoveredBy];
    end if;
    genus_list := [curves[i]`g : i in cover_keys];
    num_vals := Maximum([2*g+5 : g in genus_list]);

    abs_schofer_tab, all_cm_pts := AbsoluteValuesAtCMPoints(Xstar, curves, all_cm_pts, fs :
                                        MaxNum := num_vals, Prec := Prec,
                                        Exclude := {}, Include := Set(d_divs));
    ReduceTable(abs_schofer_tab);
    schofer_tab := ValuesAtCMPoints(abs_schofer_tab, all_cm_pts);

    if #Targets gt 0 then
        // Solve only the target covers: filter K_idxs so RationalConstraintsOnEquations
        // (which has require #ds ge 2g+3 per cover) never touches an undersolved sibling.
        schofer_tab`K_idxs := [i : i in schofer_tab`K_idxs
                                 | curves[schofer_tab`Keys_fs[i]]`W in Targets];
    end if;
    return EquationsOfCovers(schofer_tab, all_cm_pts);
end intrinsic;
```

- [ ] **Step 2: Behaviour-preserving check on the BIELLIPTIC anchor.** Write `verify_targets.m`:

```
AttachSpec("ShimuraQuotients.spec");
AttachSpec("/Users/sachihashimoto/Repos/CHIMP/CHIMP.spec");
SetVerbose("ShimuraQuotients", 0);
curves := GetHyperellipticCandidates();
D := 6; N := 29; gens := {3,29};
Xstar := rep{X : X in curves | X`D eq D and X`N eq N and IsStarCurve(X)};
W := AllALsFromGens(gens, D*N);

// Full (old) path:
cl0, ws0, ks0 := EquationsOfCovers(Xstar, curves);
k0 := rep{k : k in ks0 | curves[k]`W eq W};
f0 := HyperellipticPolynomials(cl0[Index(ks0,k0)]);

// Target-restricted path:
clT, wsT, ksT := EquationsOfCovers(Xstar, curves : Targets := {W});
kT := rep{k : k in ksT | curves[k]`W eq W};
fT := HyperellipticPolynomials(clT[Index(ksT,kT)]);

printf "targets-only returned %o cover(s) (expect 1).\n", #ksT;
printf "models equal: %o  %o\n", f0 eq fT, (f0 eq fT select "PASS" else "FAIL");
exit;
```

Run: `magma verify_targets.m 2>&1 | grep -E "PASS|FAIL|returned"`
Expected: `targets-only returned 1 cover(s)`; `models equal: true PASS`. (This is the gate before trusting any new verdict.)

- [ ] **Step 3: Commit.**

```bash
git add EquationsCovers.m verify_targets.m
git commit -m "EquationsOfCovers: add Targets to restrict num_vals + solved covers"
```

---

## Phase D — Wire drivers + launch sweeps

### Task D1: Wire `Targets` + guard + verdict rows into the three drivers

**Files:**
- Modify: `ratpts_table10.m` (`run_group`, ~line 360), `ratpts_table1.m` (`run_entry`, ~line 109), `ratpts_table6.m` (main loop, ~line 119)

- [ ] **Step 1: In each driver, before `EquationsOfCovers`, build the row's target set, call the guard, and pass `Targets`.** For `ratpts_table10.m` `run_group` (the pattern is identical for the other two — repeat it, do not factor across files):

```
// after the DIV_CUTOFF guard, before the try/EquationsOfCovers:
target_set := { AllALsFromGens(gens, D*N) : gens in gensets };
ok, req, avail := EnoughCMPointsForTargets(Xstar, curves, target_set);
if not ok then
    printf "  insufficient CM points (need %o, have %o); skipping\n", req, avail;
    for gens in gensets do
        log_row(D, N, gens, Sprintf("SKIP-insufficient-CM(need=%o,have=%o)", req, avail), "n/a");
    end for;
    return;
end if;
```
and change the cover computation to:
```
crv_list, ws, keys := EquationsOfCovers(Xstar, curves : Targets := target_set);
```
(The `Xstar` lookup already happens above the guard in each driver; reuse it. For `ratpts_table1.m` the row has a single `gens`, so `gensets := [gens]` / `target_set := {AllALsFromGens(gens,D*N)}`.)

- [ ] **Step 2: Verify the skip path logs correctly on `34,5` (Table 10 group).** Find its index, then:

Run: `magma idx:=<i_34_5> ratpts_table10.m 2>&1 | grep -E "insufficient CM"` then `grep "SKIP-insufficient-CM" ratpts_table10_results.txt | tail`
Expected: console shows `insufficient CM points (need 9, have 4)`; results file gains `SKIP-insufficient-CM(need=9,have=4)` rows for `{2,5}` and `{5,17}`, and the run finishes in well under a minute (no Borcherds compute).

- [ ] **Step 3: Verify the success path still works on `6,29` (Table 10 group) via the driver.**

Run: `magma idx:=<i_6_29> ratpts_table10.m 2>&1 | grep -E "BIELLIPTIC"`
Expected: `BIELLIPTIC` verdict, same model as the settled row in `ratpts_table10_results.txt`.

- [ ] **Step 4: Commit.**

```bash
git add ratpts_table1.m ratpts_table6.m ratpts_table10.m
git commit -m "ratpts drivers: pass Targets, add CM pre-check skip + verdict rows"
```

### Task D2: Launch the three background sweeps and report

**Files:** none (operational)

- [ ] **Step 1: Determine each table's group count.** Run each driver with no `idx` to print its group total (Table 10 prints `228 (D,N) groups`; Table 1/6 print their list lengths). Record `HI` for each.

- [ ] **Step 2: Launch sweeps in background** (timeout per group chosen from prior timings — start with `1800` s, i.e. 30 min, comfortably above the ~945 s worst non-hang):

```bash
./run_table.sh 6  1800 1 <HI6>   &   # genus-0: biggest reduction payoff, run first
./run_table.sh 1  1800 1 <HI1>   &   # genus-1
./run_table.sh 10 1800 1 <HI10>  &   # genus-2: speedup + guard skips
```

- [ ] **Step 3: Monitor and report.** Tail the `table*_run_*.log` files and summarize new verdict rows (successes, `SKIP-insufficient-CM`, `TIMEOUT`) from each `ratpts_table*_results.txt` as they land. Compare Table 6/1 rows that previously failed to confirm the reduction rescued them.

---

## Self-review notes

- **Spec coverage:** Section 1 → C1; Section 2 → B1 + D1; Section 3 → A1 + A2; Section 4 → D2; logging rows → D1 + A2. All covered.
- **Type consistency:** `Targets` / `target_set` is everywhere a `SetEnum` of `W` subgroups (`AllALsFromGens` output). `EnoughCMPointsForTargets(Xstar, curves, Targets)` signature matches its callers in D1. `cover_keys`/`target_keys` are `CurveID` indices into `curves` (consistent with `Xstar`CoveredBy` and `Keys_fs`).
- **Known caveat:** the `run_table.sh` `HI` auto-detection stub is intentionally replaced by an explicit `HI` argument (Step in A2 requires it); D2 Step 1 obtains the counts. Do not rely on the commented auto-detect.
- **Orphan processes:** if a timed-out `magma` leaves a polymake child, escalate the `perl alarm` wrapper to kill the process group (`setpgrp` + `kill -- -$pgid`); revisit only if A2/D2 show orphans.
