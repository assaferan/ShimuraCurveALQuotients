// Elliptic points of the seven open star curves, labelled by (stabiliser order, field of
// definition).  Aut(X*) permutes elliptic points preserving the order (a topological invariant)
// and, being rational (Rotger Prop 1, which survives the quotient), preserving the field.  So every
// automorphism fixes each point whose label class is a singleton, and fixes at least one point in
// each class of ODD size.  A hyperelliptic involution at genus 3 has exactly 8 fixed points, so
// nine or more odd-size classes prove X* is not hyperelliptic.
AttachSpec("ShimuraQuotients.spec"); SetColumns(0);
cs := eval Read(DATA);
for X in cs do
    DN := X`D*X`N; full := {Integers()| d : d in Divisors(DN) | GCD(d, DN div d) eq 1};
    if assigned X`IsSubhyp or X`W ne full or X`g ne 3 then continue; end if;
    printf "\n=== X_0^%o(%o)*  genus %o ===\n", X`D, X`N, X`g;
    ell := NumberOfEllipticPointsByCMOrder(X);
    tot := 0; byLabel := AssociativeArray(); byDisc := [];
    for q in Sort(SetToSequence(Keys(ell))) do
        for d in Sort(SetToSequence(Keys(ell[q]))) do
            n := ell[q][d]; if n eq 0 then continue; end if;
            tot +:= n;
            fld := "?";
            try
                K := FieldOfDefinitionOfCMPoint(X, d);
                fld := Degree(K) eq 1 select "Q" else Sprintf("deg%o:%o", Degree(K), Sprint(DefiningPolynomial(AbsoluteField(K))));
            catch e
                fld := "ERR";
            end try;
            key := Sprintf("q=%o,%o", q, fld);
            byLabel[key] := (IsDefined(byLabel, key) select byLabel[key] else 0) + n;
            Append(~byDisc, Sprintf("q=%o d=%o n=%o %o", q, d, n, fld[1..Min(28,#fld)]));
        end for;
    end for;
    print "  points:", byDisc;
    odd := 0; single := 0;
    for k in Keys(byLabel) do
        if IsOdd(byLabel[k]) then odd +:= 1; end if;
        if byLabel[k] eq 1 then single +:= 1; end if;
    end for;
    printf "  TOTAL elliptic points %o, label classes %o, singleton classes %o, ODD-SIZE classes %o%o\n",
        tot, #Keys(byLabel), single, odd, odd ge 9 select "  => NOT HYPERELLIPTIC" else "";
end for;
printf "\nELLPTS DONE\n";
quit;
