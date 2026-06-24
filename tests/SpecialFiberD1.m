// Validates the classical (D=1) special-fiber test (SpecialFiberNotHyperelliptic), including the
// composite-N reach and the Klein-four group-quotient component (push_group, |W''| up to 4):
//  (a) no contradictions -- never flags a known (geometrically) hyperelliptic X_0(N)/W;
//  (b) effective -- proves many curves, including composite N and |W''|>2 group quotients;
//  (c) the group-quotient component works (e.g. X_0(90)/<w_2,w_9>, W'' Klein-four).

procedure test_SpecialFiberD1()
    // (c) a Klein-four group-quotient component (the push_group path)
    assert SpecialFiberNotHyperelliptic(90, {Integers()|1,2,9,18});
    curves := GetHyperellipticCandidates();
    reached := 0; proven := 0; contradictions := 0; bigW := 0;
    for X in curves do
        if X`D ne 1 or X`g lt 3 then continue; end if;
        reached +:= 1;
        ok, witness := SpecialFiberNotHyperelliptic(X`N, X`W);
        if ok then
            proven +:= 1;
            if assigned X`IsHyp and X`IsHyp then
                contradictions +:= 1;
                printf "  *** CONTRADICTION: D=1 N=%o W=%o known hyperelliptic but flagged ***\n", X`N, X`W;
            end if;
            // count proofs that use a group-quotient (|W''|>2) component
            if exists{p : p in PrimeDivisors(X`N)
                      | GenusShimuraCurveQuotient(1, X`N div p, {Integers()|1}) eq 0
                        and #{w div p^Valuation(w,p) : w in X`W} gt 2} then
                bigW +:= 1;
            end if;
        end if;
    end for;
    printf "  [D1] reached=%o proven=%o (group-quotient proofs=%o) contradictions=%o\n",
           reached, proven, bigW, contradictions;
    assert contradictions eq 0;     // soundness
    assert proven gt 0;             // effectiveness
end procedure;

test_SpecialFiberD1();
