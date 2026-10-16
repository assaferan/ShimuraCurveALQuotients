// A DECISIVE CONTROL FOR THE FIXED-POINT CRITERION, which it has never had.
//
// The criterion bounds the fixed points of a hyperelliptic involution from below by the sum of
// the field degrees over the odd-size classes of elliptic points, and concludes "not
// hyperelliptic" when that sum exceeds 2g+2.  Its published control is weak: on 54 known
// hyperelliptic star curves it never fires, but nearly every sum there is 0 and the largest is 3,
// so it would pass even if the argument were wrong.
//
// EVERY GENUS-2 CURVE IS HYPERELLIPTIC.  So on a genus-2 star curve the criterion MUST return a
// sum of at most 2g+2 = 6 -- not as a matter of evidence but as a matter of fact.  A genus-2 star
// curve of large level has many elliptic points, which is exactly the regime where the criterion
// fired at X_0^6(379), so this is a demanding test and not a vacuous one.
//
// A sum above 6 here refutes the criterion (or its input data) outright.
AttachSpec("/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/ShimuraQuotients.spec");
SetColumns(0);

// the star curves of the five genus-4 members of the sixteen, plus a spread of other levels
bases := [ <2002,1>, <2030,1>, <2418,1>, <3390,1>, <6,403>,
           <1155,1>, <858,1>, <546,1>, <1122,1>, <210,13>, <6,205>, <555,1> ];
for b in bases do
    D := b[1]; N := b[2]; DN := D*N;
    if not (IsSquarefree(DN) and GCD(D,N) eq 1) then continue; end if;
    full := { Integers() | d : d in Divisors(DN) | GCD(d, DN div d) eq 1 };
    g := GenusShimuraCurveQuotient(D, N, full);
    if g ne 2 then
        printf "X_0^%o(%o)*: genus %o -- skipped (this control needs genus 2)\n", D, N, g;
        continue;
    end if;
    X := CreateShimuraQuot(D, N, full); X`g := g;
    ell := NumberOfEllipticPointsByCMOrder(X);
    classes := [* *]; bad := false;
    for q in Sort(SetToSequence(Keys(ell))) do
        for d in Sort(SetToSequence(Keys(ell[q]))) do
            n := ell[q][d]; if n eq 0 then continue; end if;
            ok := true; K := Rationals();
            try K := AbsoluteField(FieldOfDefinitionOfCMPoint(X, d)); catch e ok := false; end try;
            if not ok then bad := true; continue; end if;
            if Degree(K) ne n then
                printf "   X_0^%o(%o)*: d = %o has %o points but field degree %o\n", D,N,d,n,Degree(K);
            end if;
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
    for cl in classes do
        if IsOdd(cl[3]) then
            S +:= Degree(cl[2]);
            Append(~odd, Sprintf("(q=%o,deg=%o,n=%o,d=%o)", cl[1], Degree(cl[2]), cl[3], cl[4]));
        end if;
    end for;
    printf "X_0^%o(%o)*  genus 2, %o elliptic points in %o classes; odd classes %o\n",
        D, N, &+[ cl[3] : cl in classes ], #classes, Join(odd, " ");
    printf "    SUM = %o  against the bound 6  ==> %o%o\n", S,
        S gt 6 select "CRITERION REFUTED: it would call a genus-2 curve non-hyperelliptic"
                 else "consistent",
        bad select "   [INCOMPLETE: a field was unavailable]" else "";
end for;
printf "ELLG2 DONE\n";
quit;
