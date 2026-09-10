import "_crviso.m" : construct_crv_isomorphism;
import "_modelfile.m" : ReadModelSet;

procedure test_AllEquationsAboveCoversSingleCurve(D, N, cover_data, ws_data, curves : algebra_map := false, base_label := 0, manual_isomorphism := false, model_covers := true, model_drift_ok := false)
    // no longer needed as we now have a test for each curve
    // printf "testing equations of covers of X0*(%o;%o)...", D, N;
    assert exists(Xstar){X : X in curves | X`D eq D and X`N eq N and IsStarCurve(X)};
    covers, ws := AllEquationsAboveCovers(Xstar, curves : base_label := base_label);

    // ⚠ COUNT THE COMPARISONS ACTUALLY MADE. Without this the procedure can PASS WHILE CHECKING
    // NOTHING: the `if not is_def then continue` below silently skips every cover key that
    // AllEquationsAboveCovers did not produce, so a base whose re-derived W keys stop matching the
    // expected ones -- because a guard changed, or the pipeline started deferring a cover -- turns
    // green instead of red. That failure mode is invisible in the output, and this file is the
    // shared helper behind all 27 X0_D_N.m tests, so one blind spot here is 27 blind spots.
    // A passing check is not evidence until you know it could have failed.
    n_curve_cmp := 0;      // curve isomorphism assertions actually executed
    n_ws_cmp := 0;         // Atkin-Lehner involution assertions actually executed
    matched_Ws := {};      // cover_data keys that were actually reached

    for label in Keys(covers) do
        X := curves[label];
        is_def, datum := IsDefined(cover_data, X`W);
        if not is_def then continue; end if;
        Include(~matched_Ws, X`W);
        C_ex, scales := Explode(datum);
        P<[x]> := AmbientSpace(C_ex);
        for base in Keys(covers[label]) do
            C := covers[label][base];
            n_curve_cmp +:= 1;
            if manual_isomorphism then
                if algebra_map then
                    phi := scales;
                else
                    phi := map<C -> C_ex | Eltseq(Vector(x)*ChangeRing(scales, Universe(x)))>;
                end if;
                is_isom := IsIsomorphism(phi);
            elif (Type(C) ne CrvHyp) and (Type(C_ex) ne CrvHyp)
                 and (#DefiningPolynomials(C) eq 2) and (#DefiningPolynomials(C_ex) eq 2) then
                // ⚠ CRV PAIR: NEVER call IsIsomorphic here. Its cost tracks PRESENTATION, not
                // genus -- a genus-7 hyperelliptic curve settles in 0.06 s while the genus-3 CRV
                // pair at 14_3 runs >50 min and the genus-5 one at 26_3 >1 h (tests/IsoScreen.m).
                // That is the real reason four tests pinned a coordinate matrix, and pinning is
                // brittle: under CMNONCOPRIME=1 the pipeline re-presents 10_13's curve and the
                // hardcoded map stops being a map at all.
                // CONSTRUCT the isomorphism instead (tests/_crviso.m): take the Mobius map from
                // the hyperelliptic y-quotient, require it to carry both sides by constant
                // squares, then let IsIsomorphism certify the result. Still a PROOF -- an
                // explicit map is exhibited and checked -- and it runs in hundredths of a second.
                is_isom, phi := construct_crv_isomorphism(C, C_ex);
            else
                is_isom, phi := IsIsomorphic(C, C_ex);
            end if;
            assert is_isom;
            ws_def, ws_ex := IsDefined(ws_data, X`W);
            if not ws_def then continue; end if;

            // WHICH isomorphism, not just whether one exists.
            // phi is reused below to conjugate the Atkin-Lehner involutions, and IsIsomorphic
            // returns an ARBITRARY element of Isom(C, C_ex) -- a torsor under Aut(C_ex). Here
            // Aut is essentially the Atkin-Lehner group itself (measured: #Aut = 8 for 10_13's
            // W={1} at genus 3, 4 for 26_1's at genus 2, matching the AL group orders), so an
            // unlucky choice PERMUTES THE INVOLUTION LABELS and the check below fails on a
            // perfectly correct model. That is why four tests pin a coordinate matrix by hand.
            //
            // Pinning a matrix is brittle: under CMNONCOPRIME=1 the pipeline re-presents 10_13's
            // curve in different coordinates and the hardcoded map stops being a map at all.
            // So instead search the torsor for one that intertwines EVERY LABELLED involution
            // SIMULTANEOUSLY -- w_m must go to w_m, not to some other involution. That keeps the
            // full strength of the manual check (each named involution is still verified against
            // an explicitly exhibited isomorphism) while surviving re-presentation.
            // ⚠ PAY FOR THE SEARCH ONLY WHEN THE FIRST MAP FAILS. Computing AutomorphismGroup
            // unconditionally, once per (cover, base), is ruinous: measured at 10_13 it took the
            // test from 870 s to over 71 minutes and still climbing -- a ~5x regression on the
            // common case, where the phi that IsIsomorphic returned already works. So try phi
            // first and fall back to the torsor only if it does not intertwine the involutions.
            function ws_ok(psi)
                for Q in Keys(ws_ex) do
                    w_alg := AlgebraMap(psi)*AlgebraMap(ws[label][base][Q])*AlgebraMap(psi^(-1));
                    phi1 := map< C_ex -> C_ex | [w_alg(x[j]) : j in [1..#x]]>;
                    phi2 := map< C_ex -> C_ex | Eltseq(Vector(x)*ChangeRing(ws_ex[Q], Universe(x)))>;
                    if phi1 ne phi2 then return false; end if;
                end for;
                return true;
            end function;

            found_phi := ws_ok(phi);
            n_iso_tried := 1;
            if (not found_phi) and (not manual_isomorphism) then
                try
                    Aut, mAut := AutomorphismGroup(C_ex);
                    for a in Aut do
                        n_iso_tried +:= 1;
                        if ws_ok(phi*mAut(a)) then found_phi := true; break; end if;
                    end for;
                catch e
                    ;   // no computable automorphism group: phi was the only candidate
                end try;
            end if;
            n_ws_cmp +:= #Keys(ws_ex);
            error if not found_phi,
                Sprintf("X0^%o(%o) cover %o: no isomorphism to the expected curve intertwines all "
                        * "%o labelled Atkin-Lehner involution(s). %o candidate map(s) tried "
                        * "(Isom = Aut(C_ex) o phi). The curves ARE isomorphic -- what fails is "
                        * "that no identification matches the involution LABELLING.",
                        D, N, Sort(SetToSequence(X`W)), #Keys(ws_ex), n_iso_tried);
        end for;
    end for;

    // THE GUARD. Zero curve comparisons means nothing above was verified, however green the run
    // looks. Fail loudly and say what was expected against what was produced, so the reader can
    // tell "the pipeline stopped emitting this cover" from "the expected key is written wrong".
    // ⚠ TWO DIFFERENT CAUSES, and saying which one matters. Either no expected key matched a
    // produced cover at all, or a key DID match but `covers[label]` carried no bases, so the inner
    // loop never ran. The first version of this message reported only the first cause and
    // misdiagnosed X0_10_19 in CI, where `[1]` is both expected AND produced but has zero bases.
    error if n_curve_cmp eq 0,
        Sprintf("X0^%o(%o): NO EVIDENCE -- the test made ZERO curve comparisons, so it verified "
                * "nothing.\n  cover_data expects W in %o\n  AllEquationsAboveCovers produced W in "
                * "%o\n  %o\n  (%o expected key(s) matched a produced cover; a matched key still "
                * "yields no comparison when it has no bases -- i.e. no equation was found over "
                * "anything it covers.)",
                D, N, {Sort(SetToSequence(W)) : W in Keys(cover_data)},
                {Sort(SetToSequence(curves[l]`W)) : l in Keys(covers)},
                IsEmpty(matched_Ws)
                    select "No expected key matched a produced cover, so every comparison was skipped."
                    else "Keys MATCHED but produced no bases, so there was nothing to compare against.",
                #matched_Ws);

    // Expected covers that were never reached are NOT fatal -- a cover may legitimately be
    // deferred on a given run -- but they are silent, so say so. If this ever prints for a test
    // that is supposed to be exhaustive, that test is weaker than it looks.
    unmatched := {W : W in Keys(cover_data)} diff matched_Ws;
    if not IsEmpty(unmatched) then
        vprintf ShimuraQuotients, 1:
            "\tX0^%o(%o): %o of %o expected cover(s) were never produced, so they went unchecked: %o\n",
            D, N, #unmatched, #Keys(cover_data), {Sort(SetToSequence(W)) : W in unmatched};
    end if;
    // ⇒ THE RE-DERIVATION GAP, and why this second pass exists.
    // Everything above compares only the cover keys someone hand-wrote into cover_data, and the
    // `if not is_def then continue` at the top drops every other key IN SILENCE. Counted
    // 2026-09-09: 128 hand-written keys against 309 populated keys in data/models/ -- nine bases
    // were checking 1 of 15, eleven 1 of 4. The MODELS themselves are well checked (ModelChecks'
    // independent trace-formula point counts, ~190 Guo-Yang quotient comparisons); what was thin
    // is the claim that the PIPELINE STILL PRODUCES THEM.
    //
    // The expensive part -- AllEquationsAboveCovers -- is already paid for above, so cross-checking
    // EVERY committed cover key against `covers` costs only isomorphism tests. That is the same
    // comparison tests/_offline/ModelRegen.m makes; this reuses the run instead of paying for a
    // second one, which is what lets it be a CI check rather than an offline one.
    //
    // ⚠ IT IS A DRIFT CHECK, NOT A VALIDATION. It says "current code still produces this", not
    // "this is correct" -- the committed file is what the pipeline itself wrote. Correctness comes
    // from ModelChecks and the Guo-Yang oracles. So do NOT delete hand-written cover_data entries
    // in favour of it: those are Guo-Yang's PUBLISHED equations, which are external evidence, and
    // they are also the only entries carrying the labelled involutions (ws_data).
    n_model_cmp := 0; n_model_skip := 0; model_missing := []; model_noniso := [];
    mc_have_file := false;
    if model_covers then
        mc_have_file, mc_stored := ReadModelSet(D, N);
    end if;
    if mc_have_file then
        // Aggregate exactly as vvdata/weyl-campaign/genmodels.m does when it WRITES the file:
        // one entry per (cover label, base), keyed by Sort(W). Any other aggregation compares a
        // different object than the one the file records.
        mc_fresh := AssociativeArray();
        for mc_lab in Keys(covers) do
            mc_Wk := Sort([Integers()| w : w in curves[mc_lab]`W]);
            if not IsDefined(mc_fresh, mc_Wk) then mc_fresh[mc_Wk] := [* *]; end if;
            for mc_b in Keys(covers[mc_lab]) do
                Append(~mc_fresh[mc_Wk], covers[mc_lab][mc_b]);
            end for;
        end for;
        for mc_k in Keys(mc_stored) do
            if #mc_stored[mc_k] eq 0 then continue; end if;   // an empty entry carries no claim
            mc_key := Sort([Integers()| w : w in mc_k]);
            mc_ok, mc_list := IsDefined(mc_fresh, mc_key);
            if (not mc_ok) or (#mc_list eq 0) then
                Append(~model_missing, mc_key); continue;
            end if;
            // ⚠ MATCH AS A MULTISET -- each fresh cover absorbs at most one committed entry.
            // Without consuming matches, a key holding several curves (one per base) passes when
            // some were LOST, because two committed entries both match the same survivor. That is
            // exactly how 22_3's [1,66] hid a 3 -> 2 loss in ModelRegen's first draft.
            mc_used := {};
            for mc_e in mc_stored[mc_k] do
                if Type(mc_e[2]) eq MonStgElt then
                    // CRV entry: the file stores the defining polynomials as strings but NOT the
                    // ambient weights, so the curve cannot be rebuilt from it. Those keys are
                    // covered by the hand-written CRV data in the X0_*.m file instead.
                    n_model_skip +:= 1; continue;
                end if;
                // ⚠ A MODEL ENTRY MAY BE <genus, f, h>, meaning y^2 + h*y = f. Dropping h gives a
                // DIFFERENT curve of the same genus -- 9 entries across 7 files have one, and that
                // is what produced a false "defect" report against models_87_1 on 2026-09-08.
                if (#mc_e ge 3) and (Type(mc_e[3]) eq RngUPolElt) and (mc_e[3] ne 0) then
                    mc_Cst := HyperellipticCurve(mc_e[2], mc_e[3]);
                else
                    mc_Cst := HyperellipticCurve(mc_e[2]);
                end if;
                mc_found := false;
                for mc_i -> mc_c in mc_list do
                    if mc_i in mc_used then continue; end if;
                    if Type(mc_c) ne CrvHyp then continue; end if;
                    if Genus(mc_c) ne Genus(mc_Cst) then continue; end if;
                    if IsIsomorphic(mc_c, mc_Cst) then
                        mc_found := true; Include(~mc_used, mc_i); break;
                    end if;
                end for;
                n_model_cmp +:= 1;
                if not mc_found then Append(~model_noniso, mc_key); end if;
            end for;
        end for;

        // THE SAME GUARD AS ABOVE, for the same reason: a pass proves nothing until it could have
        // failed. Zero comparisons here means the model file was read but nothing in it was
        // comparable, which is a silent gap, not a success.
        error if n_model_cmp eq 0,
            Sprintf("X0^%o(%o): the model file was read (%o key(s)) but produced ZERO cover "
                    * "comparisons (%o CRV entr(ies) skipped) -- nothing was re-derived.",
                    D, N, #Keys(mc_stored), n_model_skip);

        // ⚠ TWO DIFFERENT FAILURES, AND ONLY ONE OF THEM IS EVER TOLERABLE.
        // A test that pins a non-zero base_label legitimately produces FEWER keys: EquationsCovers.m
        // :1061 gates EquationsByRebase on `base_label eq 0`, so the covers a default run fills by
        // rebasing the Hauptmodul stay deferred here. Measured at 10_13 (base_label 4069), whose
        // three missing keys are exactly the ones that run's own log reports as "leaving it
        // deferred". 26_3 is the case that confirms the reading rather than excusing it: its model
        // was GENERATED at base_label 8103 (data/models/PROVENANCE.md) and its test pins the same
        // value, so it must match.
        //
        // But a key the pipeline DOES produce must still be the committed curve, whatever base_label
        // was pinned. Silencing both with one flag would hide the failure that actually matters, so
        // model_drift_ok covers MISSING keys ONLY.
        error if not IsEmpty(model_noniso),
            Sprintf("X0^%o(%o): cover key(s) produced but NOT ISOMORPHIC to the curve committed in "
                    * "data/models/models_%o_%o.m: %o\n"
                    * "  (%o comparison(s) made, %o CRV entr(ies) skipped.)\n"
                    * "  This is NEVER tolerated -- model_drift_ok covers missing keys only. Note a "
                    * "re-presented but isomorphic curve does NOT reach here: the comparison is up "
                    * "to isomorphism, so this says the pipeline now builds a DIFFERENT curve.",
                    D, N, D, N, model_noniso, n_model_cmp, n_model_skip);
        error if (not model_drift_ok) and not IsEmpty(model_missing),
            Sprintf("X0^%o(%o): cover key(s) in data/models/models_%o_%o.m NOT PRODUCED AT ALL: %o\n"
                    * "  (%o comparison(s) made, %o CRV entr(ies) skipped.)\n"
                    * "  A DRIFT report, not a correctness verdict: the committed model may still be "
                    * "right (ModelChecks and the Guo-Yang oracles are what judge that). If this test "
                    * "pins a non-zero base_label, that is the expected cause and model_drift_ok is "
                    * "the answer; otherwise the pipeline has stopped producing a cover it once did.",
                    D, N, D, N, model_missing, n_model_cmp, n_model_skip);
    end if;

    vprintf ShimuraQuotients, 1:
        "\tX0^%o(%o): %o curve comparison(s), %o involution comparison(s), %o/%o expected covers "
        * "matched; %o committed model cover(s) re-derived (%o CRV skipped)\n",
        D, N, n_curve_cmp, n_ws_cmp, #matched_Ws, #Keys(cover_data), n_model_cmp, n_model_skip;
    return;
end procedure;
