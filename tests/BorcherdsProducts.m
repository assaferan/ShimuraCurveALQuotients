
procedure test_AllEquationsAboveCoversSingleCurve(D, N, cover_data, ws_data, curves : algebra_map := false, base_label := 0, manual_isomorphism := false)
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
    error if n_curve_cmp eq 0,
        Sprintf("X0^%o(%o): NO EVIDENCE -- the test made ZERO curve comparisons, so it verified "
                * "nothing.\n  cover_data expects W in %o\n  AllEquationsAboveCovers produced W in "
                * "%o\n  No expected key matched a produced cover, so every comparison was skipped.",
                D, N, {Sort(SetToSequence(W)) : W in Keys(cover_data)},
                {Sort(SetToSequence(curves[l]`W)) : l in Keys(covers)});

    // Expected covers that were never reached are NOT fatal -- a cover may legitimately be
    // deferred on a given run -- but they are silent, so say so. If this ever prints for a test
    // that is supposed to be exhaustive, that test is weaker than it looks.
    unmatched := {W : W in Keys(cover_data)} diff matched_Ws;
    if not IsEmpty(unmatched) then
        vprintf ShimuraQuotients, 1:
            "\tX0^%o(%o): %o of %o expected cover(s) were never produced, so they went unchecked: %o\n",
            D, N, #unmatched, #Keys(cover_data), {Sort(SetToSequence(W)) : W in unmatched};
    end if;
    vprintf ShimuraQuotients, 1:
        "\tX0^%o(%o): %o curve comparison(s), %o involution comparison(s), %o/%o expected covers matched\n",
        D, N, n_curve_cmp, n_ws_cmp, #matched_Ws, #Keys(cover_data);
    return;
end procedure;
