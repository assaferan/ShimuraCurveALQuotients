// tests/Genus1Classes.m
//
// THE GENUS-1 ANALOGUE OF tests/ConicClasses.m -- AND IT IS RED-BY-RECORD, NOT GREEN.
//
// WHY THIS EXISTS.  The entries stored under ONE (D, N, W) key are the SAME quotient X/W computed
// over DIFFERENT BASE covers -- `all_eqns[k][base]` in EquationsCovers.m -- so they must be
// isomorphic over Q.  tests/ConicClasses.m enforces exactly that at genus 0, by Brauer class, and
// it PASSES: 236 conics, 41 multi-entry keys, 0 disagreements.  Nothing enforced it at genus 1.
//
// Measured 2026-09-15, the first time anyone looked: of 23 multi-entry genus-1 keys, **20 hold
// entries that are not even GL2-EQUIVALENT** -- genuinely different curves filed as one.
//
// ⚠ AND THEY ALL SHARE A JACOBIAN.  Cross-checked with an invariant that does not use
// IsGL2Equivalent at all: all 14 keys sampled have entries whose Jacobians carry the SAME Cremona
// label, 14 for 14.  Inequivalent quartics with one Jacobian are DIFFERENT TORSORS of it -- so this
// is one mechanism landing on an H^1 class, not twenty independent transcription errors.  The
// suspect is the y^2 sign/scale resolution (find_y2_signs, the per-disc signs), since a choice made
// differently per base would move the torsor while leaving the Jacobian alone -- but that is an
// UNTESTED hypothesis and this repo's history is unkind to untested single-cause stories.
//
// ⚠ WHY THE VIOLATIONS ARE RECORDED RATHER THAN DELETED.  Two of them (6_5 W=[1] entry 1, 6_13
// W=[1] entry 2) are arbitrated by Gonzalez-Rotger and pinned in tests/GonzalezRotger.m's
// KNOWN_TORSOR_DRIFT.  For the other 18 there is NO oracle, so deleting an entry would be guessing
// which base got it right.  Recording them makes the class visible in CI and stops it growing.
//
// ⚠ WHAT THIS DOES NOT CATCH, same caveat as ConicClasses.m: a key whose entries are ALL wrong in
// the same way passes.  This is an internal-consistency check, not an oracle.
//
// THE TEST.  For y^2 = f with deg f <= 2g+2, two entries are the same curve over Q iff there is a
// GL2(Q) transformation [a,b,c,d] and a RATIONAL lambda with
//      f_2(x) = lambda^2 * (c x + d)^(2g+2) * f_1((a x + b)/(c x + d)).
// ⚠ The lambda^2 is not decoration: IsGL2Equivalent decides equivalence MODULO ANY SCALAR, and
// y^2 = f curves are isomorphic only when that scalar is a SQUARE.  A non-square constant is a
// different torsor -- precisely what this test is for.
// ⚠ Asymmetry, deliberate: a square constant PROVES isomorphism; finding none does NOT prove
// non-isomorphism, because IsGL2Equivalent does not promise the full orbit.  So a pass is a proof,
// and a failure is "could not be proved", which is why the violations are a recorded list and not
// an assertion that these curves differ.
//
// NB: top-level statements, not a procedure -- run_tests.m executes tests via `eval`, and an `eval`
// inside a procedure that closes over an outer variable segfaults Magma 2.29 (same trap as
// ModelChecks.m, ModelRegen.m and ConicClasses.m).  No pipeline is run: arithmetic on committed data.

g1_files := [g1_f : g1_f in Split(Pipe("ls data/models/models_*.m 2>/dev/null", ""), "\n") | #g1_f gt 0];
error if #g1_files eq 0, "Genus1Classes: no model files found under data/models/";

g1_P<g1x> := PolynomialRing(Rationals());

// The 20 pairs that cannot be proved isomorphic today: <base, W-key, entry index>.
g1_KNOWN := {
    <"10_13", [1,10], 2>, <"10_13", [1,13], 2>, <"10_13", [1,130], 2>,
    <"10_3",  [1,10], 2>, <"10_3",  [1,15], 2>, <"10_3",  [1,6],   2>,
    <"10_7",  [1,14], 2>, <"10_7",  [1,2],  2>,
    <"14_3",  [1,2],  2>, <"14_5",  [1,2],  2>, <"15_2", [1,3],   2>,
    <"22_7",  [1,154],2>, <"26_5",  [1,26], 2>,
    <"6_13",  [1],    2>, <"6_13",  [1,39], 2>, <"6_13", [1,26],  2>,
    <"6_17",  [1,2],  2>, <"6_23",  [1,138],2>, <"6_5",  [1],     2>,
    <"6_71",  [1,426],2>
};

g1_keys := 0;      // multi-entry genus->=1 keys -- the ones carrying a claim
g1_cmp  := 0;      // pairwise comparisons actually made
g1_ok   := 0;      // pairs PROVED isomorphic
g1_new  := [];     // violations not in g1_KNOWN  -- must be empty
g1_fixed := [];    // recorded violations that now PASS -- good news, but the record must be updated

for g1_fn in g1_files do
    g1_base := Substring(g1_fn, 20, #g1_fn - 21);
    g1_models := eval (Read(g1_fn) cat "\nreturn models;");
    for g1_k in Keys(g1_models) do
        // ⚠ skip CRV paired presentations: those are stored as a string, not a polynomial
        g1_es := [g1_e : g1_e in g1_models[g1_k] | Type(g1_e[2]) ne MonStgElt];
        if #g1_es lt 2 then continue; end if;             // a single entry asserts nothing here
        g1_g := g1_es[1][1];
        if g1_g lt 1 then continue; end if;               // genus 0 is ConicClasses.m's job
        error if exists{g1_e : g1_e in g1_es | g1_e[1] ne g1_g},
              Sprintf("Genus1Classes: %o W=%o holds entries of DIFFERENT genus %o -- they cannot be "
                      * "the same curve", g1_base, Sprint(g1_k), [g1_e[1] : g1_e in g1_es]);
        g1_n := 2*g1_g + 2;
        g1_fs := [];
        for g1_e in g1_es do
            // ⚠ an entry may carry an h term: y^2 + h y = f  ->  y^2 = f + h^2/4.  Dropping h gives
            // a DIFFERENT curve of the same genus (9 entries across 7 files do carry one).
            g1_ff := g1_e[2];
            if g1_e[3] ne 0 then g1_ff := g1_ff + g1_e[3]^2/4; end if;
            Append(~g1_fs, g1_ff);
        end for;
        error if exists{g1_ff : g1_ff in g1_fs | Degree(g1_ff) gt g1_n},
              Sprintf("Genus1Classes: %o W=%o claims genus %o but an entry has degree > %o",
                      g1_base, Sprint(g1_k), g1_g, g1_n);
        g1_keys +:= 1;
        for g1_i in [2..#g1_fs] do
            g1_cmp +:= 1;
            g1_good := false;
            g1_eq, g1_Ts := IsGL2Equivalent(g1_fs[1], g1_fs[g1_i], g1_n);
            if g1_eq then
                for g1_T in g1_Ts do
                    g1_a, g1_b, g1_c, g1_d := Explode(g1_T);
                    g1_den := g1_c*g1x + g1_d;
                    if g1_den eq 0 then continue; end if;
                    g1_num := g1_den^g1_n * Evaluate(g1_fs[1], (g1_a*g1x + g1_b)/g1_den);
                    if g1_num eq 0 or not IsCoercible(Rationals(), g1_fs[g1_i]/g1_num) then continue; end if;
                    g1_sq, g1_lam := IsSquare(Rationals()!(g1_fs[g1_i]/g1_num));
                    if g1_sq then
                        // certify the identity itself, not the search's say-so
                        assert g1_fs[g1_i] eq g1_lam^2 * g1_den^g1_n
                                              * Evaluate(g1_fs[1], (g1_a*g1x + g1_b)/g1_den);
                        g1_good := true; break;
                    end if;
                end for;
            end if;
            g1_tag := <g1_base, [Integers()|g1_u : g1_u in g1_k], g1_i>;
            if g1_good then
                g1_ok +:= 1;
                if g1_tag in g1_KNOWN then Append(~g1_fixed, Sprint(g1_tag)); end if;
            elif g1_tag notin g1_KNOWN then
                Append(~g1_new, Sprintf("%o W=%o entry %o (%o)", g1_base, Sprint(g1_k), g1_i,
                    g1_eq select "GL2-equivalent but the constant is a NON-SQUARE -- different torsor"
                            else "not even GL2-equivalent -- a different curve"));
            end if;
        end for;
    end for;
end for;

printf "  %o multi-entry genus->=1 key(s), %o pairwise comparison(s): %o proved isomorphic, "
       * "%o recorded as unproved\n", g1_keys, g1_cmp, g1_ok, g1_cmp - g1_ok;

// ⚠ COUNT THE COMPARISONS.  A passing check is not evidence until you know it could have failed.
error if g1_keys lt 23,
      Sprintf("Genus1Classes: only %o multi-entry key(s), expected >= 23 -- did entry parsing "
              * "regress, or did keys stop being read?", g1_keys);
error if g1_cmp lt 23,
      Sprintf("Genus1Classes: only %o comparison(s), expected >= 23 -- the check may be vacuous", g1_cmp);
error if g1_ok lt 3,
      Sprintf("Genus1Classes: only %o pair(s) PROVED isomorphic, expected >= 3 (10_7 W=[1,7], "
              * "21_2 W=[1,2], 6_13 W=[1,6]) -- the positive side has regressed", g1_ok);

error if #g1_new gt 0,
      Sprintf("Genus1Classes: %o NEW key(s) hold genus->=1 entries that are not the same curve over "
              * "Q, so one of them is wrong: %o", #g1_new, g1_new);

// ⚠ This firing is GOOD NEWS, not a regression: a recorded violation now proves isomorphic, which
// is what fixing the mechanism looks like.  Remove it from g1_KNOWN and say so in HANDOFF.md.
error if #g1_fixed gt 0,
      Sprintf("Genus1Classes: %o recorded violation(s) now PASS -- the underlying defect has been "
              * "fixed for them.  This is good news; update g1_KNOWN: %o", #g1_fixed, g1_fixed);
