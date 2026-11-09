// The quotients that the filters of the committed pipeline data leave undecided (no stage before
// the certificates proved or ruled out hyperellipticity) are exactly the 34 of data/open_34.m, the list of
// paper/hyperelliptic-candidates.tex, together with the five of its second list.  The 34 were obtained
// independently of the filters, from Noether's count on the differentials, so the test fails if a stage
// stops deciding a quotient it used to decide, or decides one of the 34.  After the certificate stage,
// exactly the listed quotients with a certificate are decided, by it.

open_34, open_gaps := eval (Read("data/open_34.m") cat "\nreturn open_34, open_gaps;");
assert #open_34 eq 34 and #open_gaps eq 5;
the34 := {<e[1], e[2], {Integers() | m : m in e[3]}> : e in open_34};
gaps := {<e[1], e[2], {Integers() | m : m in e[3]}> : e in open_gaps};
assert #the34 eq 34 and #gaps eq 5 and IsEmpty(the34 meet gaps);
expected := the34 join gaps;

// before the certificate stage: exactly the 34 and the gaps are undecided
curves := eval Read("data/curves_after_FilterByNonALInvolutions.dat");
undecided := [c : c in curves | not assigned c`IsSubhyp];
assert &and[c`g ge 3 : c in undecided];
found := {<c`D, c`N, c`W> : c in undecided};
printf "  undecided before the certificates: %o; in data/open_34.m: %o + %o\n", #undecided, #the34, #gaps;
if found ne expected then
    printf "  undecided but not in the list: %o\n", found diff expected;
    printf "  in the list but decided: %o\n", expected diff found;
end if;
assert #undecided eq #found and found eq expected;

// after it: the listed curves with a certificate are decided by it (hyperelliptic, or not for shapes N, T and V),
// nothing else moves
shape := AssociativeArray();
for k in expected do
    for c in ReadHyperellipticCertificates(k[1], k[2]) do
        if {Integers() | m : m in c["W"]} eq k[3] then shape[k] := c["shape"]; end if;
    end for;
end for;
certified := Keys(shape);
final := eval Read("data/curves_after_UpdateCurves8.dat");
left := {<c`D, c`N, c`W> : c in final | not assigned c`IsSubhyp};
printf "  certified: %o; undecided at the end: %o\n", #certified, #left;
assert left eq expected diff certified;
for c in final do
    k := <c`D, c`N, c`W>;
    if k in certified then
        assert "HyperellipticCertificate" in c`TestInWhichProved;
        hyp := shape[k] notin {"N", "T", "V"};                                    // shapes N, T, V: not hyperelliptic
        assert c`IsSubhyp eq hyp and c`IsHyp eq hyp;
    end if;
end for;

// each W is a group of Atkin-Lehner involutions of X_0^D(N): exact divisors m of DN, closed
// under w_m w_n = w_{mn/gcd(m,n)^2}
for e in expected do
    DN := e[1] * e[2];
    W := e[3];
    assert &and[DN mod m eq 0 and GCD(m, DN div m) eq 1 : m in W];
    assert &and[m * n div GCD(m, n)^2 in W : m, n in W];
end for;
