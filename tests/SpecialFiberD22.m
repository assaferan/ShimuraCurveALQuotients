// Validates the D=22 special-fiber test (SpecialFiberNotHyperellipticD22), the general CM-value
// method for a genus-0 quaternionic discriminant with no closed-form supersingular polynomial
// (star (2,2,3,3,3,3), six cone points).  Checks:
//  (a) no contradictions -- never flags a known (geometrically) hyperelliptic D=22 curve;
//  (b) effective -- proves a positive number of D=22 curves on the covered primes.
// Coverage is limited to the primes the rational CM discriminants reach (no Heun fallback), so the
// intrinsic skips uncovered primes; this test exercises the covered ones via the real curve records.

procedure test_SpecialFiberD22()
    curves := GetHyperellipticCandidates();
    // reachable W'' : the three intermediate quotients and the full star
    inter := {{Integers()|1,2},{Integers()|1,11},{Integers()|1,22},{Integers()|1,2,11,22}};
    reached := 0; proven := 0; contradictions := 0;
    for X in curves do
        if X`D ne 22 or X`g lt 3 then continue; end if;
        reachable := false;
        for p in PrimeDivisors(X`N) do
            if 22 mod p ne 0 and X`N div p eq 1 then
                Wpp := {w div p^Valuation(w,p) : w in X`W};
                if Wpp in inter then reachable := true; end if;
            end if;
        end for;
        if not reachable then continue; end if;
        reached +:= 1;
        if SpecialFiberNotHyperellipticD22(X`N, X`W) then
            proven +:= 1;
            if assigned X`IsHyp and X`IsHyp then
                contradictions +:= 1;
                printf "  *** CONTRADICTION: D=22 N=%o W=%o known hyperelliptic but flagged ***\n", X`N, X`W;
            end if;
        end if;
    end for;
    printf "  [D22] reached=%o proven=%o contradictions=%o\n", reached, proven, contradictions;
    assert contradictions eq 0;     // soundness
    assert proven gt 0;             // effectiveness (on the covered primes)
end procedure;

test_SpecialFiberD22();
