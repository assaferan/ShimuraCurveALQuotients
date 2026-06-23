// Validates the D=10 CM-value special-fiber test (SpecialFiberNotHyperellipticD10):
//  (a) no contradictions -- it never flags a known (geometrically) hyperelliptic D=10 curve as
//      non-hyperelliptic;
//  (b) it is effective -- it proves a positive number of D=10 intermediate-quotient curves
//      non-hyperelliptic, using special-fiber primes beyond the D=6 hypergeometric range.
// Curves are the genuine pipeline records (X_0(10,p)/W with W'' an intermediate quotient).

procedure test_SpecialFiberD10()
    curves := GetHyperellipticCandidates();
    inter := {{Integers()|1,2},{Integers()|1,5},{Integers()|1,10}};
    reached := 0; proven := 0; contradictions := 0; maxp := 0;
    for X in curves do
        if X`D ne 10 or X`g lt 3 then continue; end if;
        // restrict to curves the test can reach: N = p prime, W'' intermediate
        reachable := false;
        for p in PrimeDivisors(X`N) do
            if 10 mod p ne 0 and X`N div p eq 1 then
                Wpp := {w div p^Valuation(w,p) : w in X`W};
                if Wpp in inter then reachable := true; end if;
            end if;
        end for;
        if not reachable then continue; end if;
        reached +:= 1;
        ok := SpecialFiberNotHyperellipticD10(X`N, X`W);
        if ok then
            proven +:= 1;
            if X`N gt maxp then maxp := X`N; end if;     // N = p, the special-fiber prime
            if assigned X`IsHyp and X`IsHyp then
                contradictions +:= 1;
                printf "  *** CONTRADICTION: D=10 N=%o W=%o known hyperelliptic but flagged ***\n", X`N, X`W;
            end if;
        end if;
    end for;
    printf "  [D10] reached=%o proven=%o maxprime=%o contradictions=%o\n", reached, proven, maxp, contradictions;
    assert contradictions eq 0;     // soundness
    assert proven gt 0;             // effectiveness
end procedure;

test_SpecialFiberD10();
