// Sharpened criterion.  Label an elliptic point of X* by (stabiliser order q, field of definition K):
// Aut(X*) preserves q (topological) and K (its elements are rational, Rotger Prop 1, which survives
// the quotient).  Let S be a label class.  An involution iota acts on S, so #Fix(iota) and |S| have
// the same parity; if |S| is odd, iota fixes some P in S.  Being rational, iota then fixes every
// Galois conjugate of P, hence the whole orbit, of size [K:Q].  Distinct classes are disjoint, so
//        #Fix(iota)  >=  sum over odd-size classes of [K:Q].
// A hyperelliptic involution at genus 3 has exactly 8 fixed points, so a sum exceeding 8 proves
// X* is not hyperelliptic.
AttachSpec("ShimuraQuotients.spec"); SetColumns(0);
cs := eval Read(DATA);
for X in cs do
    DN := X`D*X`N; full := {Integers()| d : d in Divisors(DN) | GCD(d, DN div d) eq 1};
    if assigned X`IsSubhyp or X`W ne full or X`g ne 3 then continue; end if;
    ell := NumberOfEllipticPointsByCMOrder(X);
    classes := [* *];           // each: <q, K, total points, [discriminants]>
    bad := false;
    for q in Sort(SetToSequence(Keys(ell))) do
        for d in Sort(SetToSequence(Keys(ell[q]))) do
            n := ell[q][d]; if n eq 0 then continue; end if;
            ok := true; K := Rationals();
            try K := AbsoluteField(FieldOfDefinitionOfCMPoint(X, d)); catch e ok := false; end try;
            if not ok then bad := true; printf "X_0^%o(%o)*: no field at d = %o\n", X`D, X`N, d; continue; end if;
            if Degree(K) ne n then printf "X_0^%o(%o)*: WARNING d = %o has %o points but field degree %o\n", X`D, X`N, d, n, Degree(K); end if;
            placed := false;
            for i in [1..#classes] do
                if classes[i][1] eq q and IsIsomorphic(classes[i][2], K) then
                    classes[i][3] +:= n; Append(~classes[i][4], d); placed := true; break;
                end if;
            end for;
            if not placed then Append(~classes, [* q, K, n, [d] *]); end if;
        end for;
    end for;
    S := 0; odd := [];
    for c in classes do
        if IsOdd(c[3]) then S +:= Degree(c[2]); Append(~odd, Sprintf("(q=%o,deg %o,n=%o,d=%o)", c[1], Degree(c[2]), c[3], c[4])); end if;
    end for;
    printf "X_0^%o(%o)*: %o classes; odd-size classes %o; SUM OF DEGREES = %o%o%o\n",
        X`D, X`N, #classes, Join(odd, " "), S,
        S gt 8 select "  ==> NOT HYPERELLIPTIC" else "  (no conclusion, need > 8)",
        bad select "  [INCOMPLETE: a field was unavailable]" else "";
end for;
printf "ELL2 DONE\n";
quit;
