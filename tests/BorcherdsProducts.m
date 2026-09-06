
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
            for Q in Keys(ws_ex) do
                w_alg := AlgebraMap(phi)*AlgebraMap(ws[label][base][Q])*AlgebraMap(phi^(-1));
                phi1 := map< C_ex -> C_ex | [w_alg(x[j]) : j in [1..#x]]>;
                phi2 := map< C_ex -> C_ex | Eltseq(Vector(x)*ChangeRing(ws_ex[Q], Universe(x)))>;
                n_ws_cmp +:= 1;
                assert phi1 eq phi2;
            end for;
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
