// Does every committed model for a base with PUBLISHED data actually have an oracle?
//
// ⚠ WHY THIS EXISTS. Nothing in this repo asked that question until 2026-09-22, and the cost was
// concrete: `models_69_1.m` landed 2026-09-14 (`4bfb859`, the ScaleForSchofer w_1 fix) and sat for
// eight days validated by NOTHING external -- no entry in tests/GuoYangEquations.m, no X0_69_1.m --
// although Guo-Yang PUBLISH its equation. Only ModelChecks covered it, and ModelChecks is
// structural: genus self-consistency, the genus formula, Weil divisibility, point counts. Every one
// of those passes on a wrong curve of the right genus.
//
// The cause was not an oversight about 69_1. It is that the oracles are wired up by SWEEPS -- the
// X0_* test batches of 2025-11-19 and 2026-09-06/07, and whenever someone last extended the
// equation table -- and a model that arrives AFTER a sweep is never swept. That will recur every
// time a fix builds a new base, which is exactly when attention is on the fix. So the coverage
// property is checked here instead of being left to whoever remembers.
//
// WHAT THIS IS NOT. It does not check that a model is CORRECT -- the oracles do that. It checks
// that an oracle EXISTS. A base passing here is a base some other test is looking at; a base
// failing here is one nobody is.
//
// ✅ NEGATIVE-CONTROLLED 2026-09-22. Three ways it must go red, all exercised; re-run them after
// any change here, because a coverage test that cannot fail is worse than none -- it converts an
// unknown into a false assurance.
//   1. a base loses its only oracle:  mv tests/X0_94_1.m /tmp && <run> && mv back
//        -> "1 Guo-Yang base(s) have a COMMITTED MODEL and NO ORACLE: [ 94_1 ]"
//        This is the 69_1 bug exactly, replayed on a base that still has one.
//   2. an exemption goes stale:       printf '//' > tests/X0_15_4.m && <run> && rm
//        -> "[ 15_4 ] is listed as EXEMPT but now HAS an oracle"
//   3. the search pattern rots:       break the literal in eq_pattern_hit
//        -> the 51_1 canary fires.
// ⚠ CONTROL 3 IS WHY THE CANARY EXISTS. The first version of this file FAILED that control: with
// the pattern fully broken it still printed "ok", because 93_1's hardcoded special case held the
// count at 1 and the `eq 0` guard could not fire. The bug was in the guard, not the logic, and
// only the control found it.
//
// SCOPE: the repo's TWO published sources, swept separately because they are published at
// different granularities and the difference is exactly where the last gap hid.
//   PART A  Guo-Yang, 43 rows, granularity = BASE (D,N)          -- full equations
//   PART B  Gonzalez-Rotger, 31 rows, granularity = BASE OR KEY  -- genus-one curves AND quotients
// The two are COMPLETELY DISJOINT as base sets: measured 2026-09-23, GR's genus-one bases share not
// one (D,N) with Guo-Yang's 43. So neither is redundant cover for the other.
//
// ⚠ PART B EXISTS BECAUSE THE BASE WAS THE WRONG UNIT. Until 2026-09-23 this file swept bases only,
// and read "GR" as its eleven genus-one curves (Table 1). GR also publish Table 2: seventeen
// genus-one ATKIN-LEHNER QUOTIENTS X_D^(m), which are our W=[1,m] KEYS, not bases. Seven had a
// committed model and none was checked against it -- and every one of those seven sits on a base
// (39_1 55_1 62_1 69_1 77_1 94_1 178_1) that this file already reported as COVERED, because its
// W=[1] curve has a Guo-Yang equation. ⇒ A GREEN VERDICT ON A BASE SAID NOTHING ABOUT ITS QUOTIENT
// KEYS, which is the 69_1 lesson one level down: a sweep is blind to objects it does not enumerate.
// All seven now verify by an exhibited Q-isomorphism (tests/GonzalezRotger.m PART 4).

printf "Checking oracle coverage of committed models...";

// THE 43 GUO-YANG EQUATION BASES.
// ⚠ PROVENANCE, because a wrong list here silently weakens every verdict below. Determined from
// the paper source by counting EQUATION CELLS, not labels: `\multirow{1}{*}{\text}` occurs exactly
// 43 times, and each cell's curve label sits exactly two lines above it, giving 43 cells to 43
// distinct labels with no repeats. Counting labels instead returns 41 -- two rows write the label
// without braces round D (`$X^6_0(17)$`, `$X^6_0(29)$`) -- and a bare grep for a base number hits
// coefficient digits. Cross-checked against the journal (Compositio 153 (2017)), Table A.1.
GY := [ <10,11>, <10,13>, <10,19>, <10,23>, <111,1>, <119,1>, <134,1>, <146,1>, <14,3>, <14,5>,
        <159,1>, <15,2>, <15,4>, <194,1>, <206,1>, <21,2>, <22,3>, <22,5>, <26,1>, <26,3>,
        <35,1>, <38,1>, <39,1>, <39,2>, <51,1>, <55,1>, <57,1>, <58,1>, <62,1>, <69,1>,
        <6,11>, <6,17>, <6,19>, <6,29>, <6,31>, <6,37>, <74,1>, <82,1>, <86,1>, <87,1>,
        <93,1>, <94,1>, <95,1> ];
error if #GY ne 43,
    Sprintf("OracleCoverage: the Guo-Yang base list has %o entries, not the 43 equation cells "
            * "counted in the source -- fix the list before trusting any verdict below", #GY);

// DOCUMENTED EXEMPTIONS: a base that legitimately has no Guo-Yang oracle, with the reason.
// ⚠ Keep this list minimal and justified. An exemption is a hole someone decided to accept; an
// undocumented one is a hole nobody noticed.
//   15_4 -- the JOURNAL's Remark 39 states X_0^15(4) is OUTSIDE Guo-Yang's method (the normalizer
//           of the Eichler order strictly contains the Atkin-Lehner group), which is why the
//           honest denominator for reproduction is 42 and not 43. It is checked instead by
//           tests/CRV_15_4.m, which is not in the X0_ naming form.
EXEMPT := [ <15,4> ];

gy_src := Read("tests/GuoYangEquations.m");

// Covered by the equation table? The case tuples are written `<D, N, [Integers()|...`, in gy_cases,
// in gy_pairs and in the PENDING block alike, so one pattern finds all three.
// ⚠ DELIBERATELY CONSERVATIVE: this is a literal substring search, so reformatting that file makes
// a base read as UNCOVERED and turns this test RED. That is the safe direction -- a false alarm
// costs a minute, a false all-clear is what this file exists to prevent. 93_1 is hardcoded because
// it is handled in its own block rather than as a case tuple.
// ⚠ THE HARDCODED 93_1 IS KEPT OUT OF THE PATTERN COUNT, and that detail is load-bearing. The
// first version folded it in, so when a control broke the search pattern completely the count fell
// to 1 rather than 0, the `oc_eq eq 0` guard below could not fire, and the test PASSED while
// detecting nothing -- a guard incapable of failing, which is the defect this file exists to
// prevent. Pattern hits are now counted separately from the special case.
function eq_pattern_hit(D, N)
    return Index(gy_src, Sprintf("<%o, %o, [Integers()|", D, N)) gt 0;
end function;

function covered_by_equations(D, N)
    return eq_pattern_hit(D, N) or (<D,N> eq <93,1>);
end function;

// CANARY: prove the search works against a KNOWN entry before trusting any verdict it produces.
// 51_1 has been a case tuple since the file was written. If this fires, tests/GuoYangEquations.m
// was reformatted and every "covered by the equation table" verdict below is worthless.
error if not eq_pattern_hit(51, 1),
    "OracleCoverage: the equation-table search no longer matches the known entry 51_1, so it can "
    * "no longer tell coverage from absence. tests/GuoYangEquations.m was probably reformatted -- "
    * "update the pattern in eq_pattern_hit before trusting this test again.";

function covered_by_x0(D, N)
    return FileExists(Sprintf("tests/X0_%o_%o.m", D, N))
        or FileExists(Sprintf("tests/_offline/X0_%o_%o.m", D, N));
end function;

// Which bases have a committed model? Same enumeration idiom as ModelChecks.m.
oc_files := Split(Pipe("ls data/models/models_*.m 2>/dev/null", ""), "\n");
oc_files := [f : f in oc_files | #f gt 0];
error if #oc_files eq 0, "OracleCoverage: no model files found under data/models/";
have_model := {};
for f in oc_files do
    parts := Split(f, "/");
    name  := parts[#parts];
    core  := name[8..#name-2];                    // strip "models_" and ".m"
    dn    := Split(core, "_");
    if #dn ne 2 then continue; end if;
    Include(~have_model, <StringToInteger(dn[1]), StringToInteger(dn[2])>);
end for;

oc_checked := 0;      // GY bases with a model that were actually examined
oc_eq      := 0;
oc_eq_pat  := 0;      // equation-table hits found by the PATTERN alone, excluding 93_1
oc_x0      := 0;
oc_naked   := [];     // model, no oracle, not exempt  -> FAILURE
oc_nomodel := [];     // published base we cannot yet check, because nothing is built
oc_stale   := [];     // exempt but now covered -> the exemption should be retired

for b in GY do
    D, N := Explode(b);
    if <D,N> notin have_model then
        Append(~oc_nomodel, Sprintf("%o_%o", D, N));
        continue;
    end if;
    oc_checked +:= 1;
    eq_ok := covered_by_equations(D, N);
    x0_ok := covered_by_x0(D, N);
    if eq_ok then oc_eq +:= 1; end if;
    if eq_pattern_hit(D, N) then oc_eq_pat +:= 1; end if;
    if x0_ok then oc_x0 +:= 1; end if;
    if <D,N> in EXEMPT then
        if eq_ok or x0_ok then Append(~oc_stale, Sprintf("%o_%o", D, N)); end if;
        continue;
    end if;
    if not (eq_ok or x0_ok) then Append(~oc_naked, Sprintf("%o_%o", D, N)); end if;
end for;

// NON-VACUITY. If the list, the enumeration or the pattern ever stops matching reality this test
// would sail through having examined nothing -- the exact failure mode it is meant to catch.
error if oc_checked eq 0,
    "OracleCoverage: NO EVIDENCE -- zero Guo-Yang bases with a committed model were examined, so "
    * "this test verified nothing. Either data/models/ is empty or the GY list is wrong.";
error if oc_eq_pat eq 0 or oc_x0 eq 0,
    Sprintf("OracleCoverage: the coverage probes look broken -- %o base(s) matched the equation "
            * "table BY PATTERN (93_1's special case excluded on purpose, so this can reach 0) and "
            * "%o matched an X0_ test. Both should be well above zero; a literal-search pattern "
            * "that stopped matching reads as total absence of coverage.", oc_eq_pat, oc_x0);

error if #oc_stale ne 0,
    Sprintf("OracleCoverage: %o is listed as EXEMPT but now HAS an oracle. Retire the exemption -- "
            * "a stale exemption hides the next real gap.", oc_stale);

error if #oc_naked ne 0,
    Sprintf("OracleCoverage: %o Guo-Yang base(s) have a COMMITTED MODEL and NO ORACLE: %o.\n"
            * "  Guo-Yang publish an equation for each, so the model is checked only structurally "
            * "(ModelChecks), which passes on a wrong curve of the right genus.\n"
            * "  Fix by adding the published equation to tests/GuoYangEquations.m (cheap, "
            * "milliseconds) and/or a tests/X0_D_N.m re-derivation test. If the base genuinely "
            * "cannot have one, add it to EXEMPT above WITH THE REASON.",
            #oc_naked, oc_naked);

// =============================================================================================
// PART B: GONZALEZ-ROTGER. Unit is a published OBJECT -- a base for Table 1, a quotient KEY for
// Table 2 and footnote 2 -- because that is what the paper publishes and what a model can be wrong
// about independently.
//
// ⚠ THE LISTS ARE THE PAPER'S OWN, NOT A SELECTION OF WHAT WE HAPPEN TO CHECK. Lemma 3.1 states
// the 11 genus-one X_0(D,N), Lemma 4.1 the 17 genus-one non-elliptic X_D^(m), footnote 2 (p.11) the
// 3 remaining X_D^(m) that ARE elliptic over Q. Transcribing "what GonzalezRotger.m happens to
// contain" would make this test a tautology -- it would sweep its own answer.
GR_T1 := [ <14,1>, <15,1>, <21,1>, <33,1>, <34,1>, <46,1>, <6,5>, <6,7>, <6,13>, <10,3>, <10,7> ];
GR_T2 := [ <39,13>, <55,5>, <62,2>, <69,3>, <77,11>, <85,17>, <94,2>, <178,89>,
           <210,30>, <210,42>, <210,70>, <210,105>,
           <330,3>, <330,22>, <330,33>, <330,165>, <462,154> ];
GR_FN := [ <35,7>, <51,3>, <115,23> ];
error if #GR_T1 ne 11 or #GR_T2 ne 17 or #GR_FN ne 3,
    Sprintf("OracleCoverage: the Gonzalez-Rotger lists have %o/%o/%o entries, not the published "
            * "11/17/3 (Lemma 3.1, Lemma 4.1, footnote 2) -- fix them before trusting any verdict",
            #GR_T1, #GR_T2, #GR_FN);

gr_src := Read("tests/GonzalezRotger.m");

// Same DELIBERATELY CONSERVATIVE literal search as PART A: reformatting tests/GonzalezRotger.m
// makes an object read as UNCOVERED and turns this test red, which is the safe direction.
// ⚠ The "<D,m," spelling must stay unambiguous. MEASURED over all 31 rows, not eyeballed: no
// pattern is a substring of another -- e.g. "<330,3," is not inside "<330,33,", because the
// character after "<330,3" there is "3", not ",".
// ⚠ KNOWN WEAKNESS, stated rather than papered over: the search finds a pair ANYWHERE in the file,
// and "<6,5," and "<6,13," each occur three times -- in the Table 1 list, the quotient list, and
// NO_EXHIBITED_ISO. So deleting just those two Table 1 rows would still read as covered here. It is
// not worth more machinery: that deletion makes tests/GonzalezRotger.m itself go red on its own
// comparison counts, so the case is caught, just by the other test.
function gr_hit(D, m)
    return Index(gr_src, Sprintf("<%o,%o,", D, m)) gt 0;
end function;

// CANARY, for the reason PART A has one: prove the search matches a KNOWN row of each of the three
// tables before believing any absence it reports. One canary per table, because the three are
// written in different blocks and could be reformatted independently.
error if not (gr_hit(14,1) and gr_hit(39,13) and gr_hit(35,7)),
    Sprintf("OracleCoverage: the Gonzalez-Rotger search no longer matches known rows (Table 1 "
            * "(14,1): %o, Table 2 (39,13): %o, footnote (35,7): %o), so it can no longer tell "
            * "coverage from absence -- tests/GonzalezRotger.m was probably reformatted. Update "
            * "gr_hit before trusting this test again.", gr_hit(14,1), gr_hit(39,13), gr_hit(35,7));

// Does a committed model actually hold this object? For a Table 1 row that is the W=[1] key of
// (D,N); for a quotient row, the W=[1,m] key of (D,1). An absent or empty key means there is
// nothing yet to check, not a gap.
function gr_have(D, N, key)
    fn := Sprintf("data/models/models_%o_%o.m", D, N);
    if not FileExists(fn) then return false; end if;
    models := eval (Read(fn) cat "\nreturn models;");
    return IsDefined(models, key) and #models[key] gt 0;
end function;

gr_checked := 0; gr_naked := []; gr_nomodel := 0;
for t in GR_T1 do
    D, N := Explode(t);
    if not gr_have(D, N, [Integers()|1]) then gr_nomodel +:= 1; continue; end if;
    gr_checked +:= 1;
    if not gr_hit(D, N) then Append(~gr_naked, Sprintf("%o_%o W=[1]", D, N)); end if;
end for;
for t in GR_T2 cat GR_FN do
    D, m := Explode(t);
    if not gr_have(D, 1, Sort([Integers()|1, m])) then gr_nomodel +:= 1; continue; end if;
    gr_checked +:= 1;
    if not gr_hit(D, m) then Append(~gr_naked, Sprintf("%o_1 W=[1,%o]", D, m)); end if;
end for;

error if gr_checked eq 0,
    "OracleCoverage: NO EVIDENCE for Gonzalez-Rotger -- zero published GR objects with a committed "
    * "model were examined, so PART B verified nothing.";
error if #gr_naked ne 0,
    Sprintf("OracleCoverage: %o Gonzalez-Rotger object(s) have a COMMITTED MODEL and NO ORACLE: "
            * "%o.\n  Gonzalez-Rotger publish an equation (Table 1 / Table 2) or a Cremona label "
            * "(footnote 2) for each, so the model is currently checked only structurally, which "
            * "passes on a wrong curve of the right genus.\n  Fix by adding the published row to "
            * "tests/GonzalezRotger.m -- it is a transcription, and the file's self-check will "
            * "verify the row against the paper's own stated Jacobian before comparing ours.",
            #gr_naked, gr_naked);

printf " ok (Guo-Yang: %o of %o published base(s) have a model; %o via the equation table, %o via "
       * "an X0_ test; %o exempt; %o not yet built: %o.  Gonzalez-Rotger: %o of %o published "
       * "object(s) have a model, all with an oracle; %o not yet built)\n",
       oc_checked, #GY, oc_eq, oc_x0, #EXEMPT, #oc_nomodel, oc_nomodel,
       gr_checked, #GR_T1 + #GR_T2 + #GR_FN, gr_nomodel;
