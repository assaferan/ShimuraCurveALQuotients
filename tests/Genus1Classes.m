// tests/Genus1Classes.m
//
// THE GENUS-1 ANALOGUE OF tests/ConicClasses.m: entries filed under one key must describe the
// same curve, so their JACOBIANS must agree.
//
// WHY THIS EXISTS.  The entries stored under ONE (D, N, W) key are the SAME quotient X/W computed
// over DIFFERENT BASE covers -- `all_eqns[k][base]` in EquationsCovers.m -- so they must be
// isomorphic over Q.  ConicClasses.m enforces that at genus 0 by Brauer class and passes (236
// conics, 41 multi-entry keys, 0 disagreements).  Nothing enforced anything at genus 1.
//
// ⚠⚠ WHAT THIS TEST MUST NOT DO, AND WHY -- READ BEFORE "STRENGTHENING" IT.
//
// The obvious check is "the two quartics must be GL2-equivalent".  THAT IS THE WRONG CRITERION, and
// it was written, committed and retracted on 2026-09-15 within the hour.  Measured then: 20 of 23
// multi-entry genus-1 keys hold entries that are not GL2-equivalent, and every one of those pairs
// has the SAME Jacobian.  That looks damning and is in fact EXPECTED:
//
//   * a degree-2 map from a genus-1 curve C to P^1 is a Q-rational degree-2 divisor class, and
//     those form a torsor under E(Q) = Pic^0(C);
//   * so when E(Q) is nontrivial, ONE CURVE carries SEVERAL INEQUIVALENT QUARTIC MODELS --
//     inequivalent as binary quartics, identical as curves;
//   * and E(Q) is nontrivial for every Jacobian in play here.  Measured:
//         30a6 [2,2]  78a2 [2,2]  30a2 [2,6]  30a3 [2]  102b3 [2,2]  130b2 [2,2]
//         154a1 rank 1 [2]   138a1 rank 1 [2]   426b1 rank 1 [2]
//     Not one trivial group.
//
// Two bases picking different degree-2 classes therefore produce inequivalent quartics of the same
// curve, which is exactly the observed signature.  Corroborating: the entries of all eight keys
// probed have IDENTICAL everywhere-local solubility profiles (real place and every prime to 47).
// ⇒ "not GL2-equivalent" is NOT evidence of a wrong model.  Do not assert on it.
//
// ⚠ AND THE SAME CORRECTION APPLIES TO tests/GonzalezRotger.m: an entry there that cannot be proved
// isomorphic to the published quartic is not thereby shown to be the wrong curve.
//
// WHAT IS CHECKED HERE, then: the Jacobian, which IS an invariant of the curve and is blind to the
// choice of degree-2 class.  ⚠ It is NECESSARY and NOT SUFFICIENT -- inequivalent torsors of one E
// share it -- exactly as ConicClasses.m is an internal-consistency check and not an oracle.  At
// genus 0 the Brauer class happens to be a COMPLETE invariant of a conic; at genus 1 the Jacobian
// is not, and that gap is real and unclosed.  A key whose entries are all wrong the same way passes.
//
// NB: top-level statements, not a procedure -- run_tests.m executes tests via `eval`, and an `eval`
// inside a procedure that closes over an outer variable segfaults Magma 2.29 (same trap as
// ModelChecks.m, ModelRegen.m and ConicClasses.m).  No pipeline is run: arithmetic on committed data.

g1_files := [g1_f : g1_f in Split(Pipe("ls data/models/models_*.m 2>/dev/null", ""), "\n") | #g1_f gt 0];
error if #g1_files eq 0, "Genus1Classes: no model files found under data/models/";

g1_P<g1x> := PolynomialRing(Rationals());

// Jacobian of y^2 = f(x), deg f <= 4, via the classical invariants (Gonzalez-Rotger Section 2, p.3).
// ⚠ NOT via Jacobian()/EllipticCurve(): those want a rational point, and these curves have none by
// construction -- that route returns ERR on both sides and would compare nothing.
function g1_jac(f)
    g1_c := [Coefficient(f, g1_i) : g1_i in [0..4]];
    g1_I := 12*g1_c[5]*g1_c[1] - 3*g1_c[4]*g1_c[2] + g1_c[3]^2;
    g1_J := 72*g1_c[5]*g1_c[3]*g1_c[1] + 9*g1_c[4]*g1_c[3]*g1_c[2] - 27*g1_c[5]*g1_c[2]^2
            - 27*g1_c[4]^2*g1_c[1] - 2*g1_c[3]^3;
    return MinimalModel(EllipticCurve([0, 0, 0, -27*g1_I, -27*g1_J]));
end function;

g1_keys := 0;      // multi-entry keys of genus >= 1 -- the ones carrying a claim
g1_cmp  := 0;      // pairwise Jacobian comparisons actually made
g1_bad  := [];

for g1_fn in g1_files do
    g1_base := Substring(g1_fn, 20, #g1_fn - 21);
    g1_models := eval (Read(g1_fn) cat "\nreturn models;");
    for g1_k in Keys(g1_models) do
        // ⚠ skip CRV paired presentations: stored as a string, not a polynomial
        g1_es := [g1_e : g1_e in g1_models[g1_k] | Type(g1_e[2]) ne MonStgElt];
        if #g1_es lt 2 then continue; end if;             // a single entry asserts nothing here
        g1_g := g1_es[1][1];
        if g1_g lt 1 then continue; end if;               // genus 0 is ConicClasses.m's job
        error if exists{g1_e : g1_e in g1_es | g1_e[1] ne g1_g},
              Sprintf("Genus1Classes: %o W=%o holds entries of DIFFERENT genus %o -- they cannot be "
                      * "the same curve", g1_base, Sprint(g1_k), [g1_e[1] : g1_e in g1_es]);
        // Every multi-entry key of genus >= 1 in the committed data is genus 1.  If a genus-2 one
        // appears, g1_jac is the wrong invariant for it and this must stop rather than mis-compare.
        error if g1_g ne 1,
              Sprintf("Genus1Classes: %o W=%o is a multi-entry key of genus %o; only genus 1 was "
                      * "ever present, and the Jacobian invariant used here is genus-1 only",
                      g1_base, Sprint(g1_k), g1_g);
        g1_fs := [];
        for g1_e in g1_es do
            // ⚠ an entry may carry an h term: y^2 + h y = f  ->  y^2 = f + h^2/4.  Dropping h gives
            // a DIFFERENT curve of the same genus (9 entries across 7 files do carry one).
            g1_ff := g1_e[2];
            if g1_e[3] ne 0 then g1_ff := g1_ff + g1_e[3]^2/4; end if;
            error if Degree(g1_ff) gt 4,
                  Sprintf("Genus1Classes: %o W=%o claims genus 1 but an entry has degree %o",
                          g1_base, Sprint(g1_k), Degree(g1_ff));
            Append(~g1_fs, g1_ff);
        end for;
        g1_keys +:= 1;
        g1_E1 := g1_jac(g1_fs[1]);
        for g1_i in [2..#g1_fs] do
            g1_cmp +:= 1;
            g1_Ei := g1_jac(g1_fs[g1_i]);
            if not IsIsomorphic(g1_E1, g1_Ei) then
                Append(~g1_bad, Sprintf("%o W=%o entry 1 vs %o: Jacobians %o and %o",
                       g1_base, Sprint(g1_k), g1_i, aInvariants(g1_E1), aInvariants(g1_Ei)));
            end if;
        end for;
    end for;
end for;

printf "  %o multi-entry genus-1 key(s), %o pairwise Jacobian comparison(s)\n", g1_keys, g1_cmp;

// ⚠ COUNT THE COMPARISONS.  A passing check is not evidence until you know it could have failed.
error if g1_keys lt 23,
      Sprintf("Genus1Classes: only %o multi-entry key(s), expected >= 23 -- did entry parsing "
              * "regress, or did keys stop being read?", g1_keys);
error if g1_cmp lt 23,
      Sprintf("Genus1Classes: only %o comparison(s), expected >= 23 -- the check may be vacuous", g1_cmp);

error if #g1_bad gt 0,
      Sprintf("Genus1Classes: %o key(s) hold genus-1 entries with DIFFERENT Jacobians, so they are "
              * "not the same curve and one of them is wrong: %o", #g1_bad, g1_bad);

printf "  all %o multi-entry key(s) agree on the Jacobian\n", g1_keys;
