// tests/CMFieldsRegressionPin.m
//
// REGRESSION PIN -- pinned against the code's output as of 2026-09-24, consistent across
// versions 081d1ae, fddce4e and current; guarantees UNCHANGED, not CORRECT; for source-backed
// checks see tests/CMPoints.m and tests/CMFieldsOfDefinition.m.
//
// What it checks.  data/cm_fields_pin.m lists 3252 CM points (curve, d): 26 quotients
// X_0(D,N)/W with D*N squarefree (star quotients and proper subquotients, D = 1 and D > 1,
// N = 1 and N > 1) times the 151 discriminants of class number <= 4.  1133 carry a CM point and
// pin its field(s) of definition; 2119 pin "no CM point by this order" (the existence test,
// including the GCD(D, f) = 1 correction).  Each entry was kept only where the slow and fast
// field-of-definition functions agreed in all three versions and the degree function matched;
// the 674 keys where anything disagreed are listed in PIN_EXCLUDED and NOT checked here (see the
// data file's header for what disagreed).  For every pinned entry this test asserts that the
// current FieldsOfDefinitionOfCMPointFast returns the same fields UP TO ISOMORPHISM (as sets) and
// that DegreeOfFieldOfDefinitionOfCMPoint returns the pinned degree.
//
// A FAILURE HERE MEANS THE ANSWER CHANGED, not that it is wrong.  If the change is intended,
// regenerate with tools/regen-cm-fields-pin.m, review the diff line by line, and re-pin.
// Runtime ~45 s.

procedure test_CMFieldsRegressionPin(PIN, PIN_CURVES, PIN_DISCS, PIN_EXCLUDED)
    printf "Testing CM fields of definition against the regression pin...";

    // Guard on the table itself: every key of PIN_CURVES x PIN_DISCS is pinned or excluded,
    // exactly once, so a truncated or hand-edited table cannot silently shrink the check.
    pinned   := [<e[1], e[2], e[3], e[4]> : e in PIN];
    excluded := [<e[1], e[2], e[3], e[4]> : e in PIN_EXCLUDED];
    assert #Set(pinned) eq #pinned and #Set(excluded) eq #excluded;
    assert #(Set(pinned) meet Set(excluded)) eq 0;
    assert Set(pinned) join Set(excluded) eq
        {<c[1], c[2], c[3], d> : c in PIN_CURVES, d in PIN_DISCS};
    assert #PIN eq 3252 and #PIN_EXCLUDED eq 674;

    R<x> := PolynomialRing(Rationals());
    toabs := func<F | Type(F) eq FldRat select F else AbsoluteField(F)>;
    iso := func<F, G | AbsoluteDegree(F) eq AbsoluteDegree(G) and
                       (AbsoluteDegree(F) eq 1 or IsIsomorphic(F, G))>;
    Xs := AssociativeArray();
    failures := [];
    checked := 0;
    for e in PIN do
        D, N, W, d, flds, deg := Explode(e);
        if not IsDefined(Xs, <D, N, W>) then
            Xs[<D, N, W>] := CreateShimuraQuot(D, N, Set(W));
        end if;
        X := Xs[<D, N, W>];
        pinF := [* toabs(#c eq 2 select Rationals() else NumberField(R ! c)) : c in flds *];
        got := [* toabs(F) : F in FieldsOfDefinitionOfCMPointFast(X, d) *];
        // equal as sets up to isomorphism: each side's every field is isomorphic to one on the other
        same := true;
        for pair in [<pinF, got>, <got, pinF>] do
            for F in pair[1] do
                if not exists{i : i in [1..#pair[2]] | iso(F, pair[2][i])} then same := false; end if;
            end for;
        end for;
        got_deg := DegreeOfFieldOfDefinitionOfCMPoint(X, d);
        if not same or got_deg ne deg then
            Append(~failures, Sprintf("(%o, %o, %o) d = %o: pinned %o of degree %o, now %o of degree %o",
                D, N, W, d, flds, deg, [Type(F) eq FldRat select [0, 1] else Coefficients(Polredabs(DefiningPolynomial(F))) : F in got], got_deg));
        end if;
        checked +:= 1;
    end for;
    assert checked eq #PIN;
    error if #failures gt 0,
        Sprintf("%o of %o pinned CM points changed:\n%o", #failures, #PIN, Join(failures, "\n"));
    printf "Done! (%o points checked, %o with a field)\n", checked, #[e : e in PIN | #e[5] gt 0];
end procedure;

// The table is read HERE, at top level, and passed in: Magma 2.29-4 SEGFAULTS on an `eval` of this
// file made inside a procedure or function body when the test itself runs under run_tests.m's eval.
cmpin_src := Read("data/cm_fields_pin.m");
test_CMFieldsRegressionPin(eval (cmpin_src cat "\nreturn PIN;"),
                           eval (cmpin_src cat "\nreturn PIN_CURVES;"),
                           eval (cmpin_src cat "\nreturn PIN_DISCS;"),
                           eval (cmpin_src cat "\nreturn PIN_EXCLUDED;"));
