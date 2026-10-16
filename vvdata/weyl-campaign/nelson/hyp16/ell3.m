// The repository's own fixed-point criterion for star curves, applied to the three members of
// the sixteen that ARE star curves (W = the full Atkin-Lehner group).  An involution of a star
// curve is rational (Rotger, Prop. 1), so it fixes the whole Galois orbit of any point it fixes;
// an odd-size class of elliptic points labelled by (stabiliser order, field of definition) must
// therefore contain a fixed point, and the degrees of those fields add up to a LOWER BOUND for
// the number of fixed points.  A hyperelliptic involution at genus 3 has exactly 8.
AttachSpec("/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/ShimuraQuotients.spec");
SetColumns(0);
cases := [ <46,31>, <6,379>, <6,649> ];
for c in cases do
    D := c[1]; N := c[2]; DN := D*N;
    full := { Integers() | d : d in Divisors(DN) | GCD(d, DN div d) eq 1 };
    X := CreateShimuraQuot(D, N, full);
    X`g := GenusShimuraCurveQuotient(D, N, full);
    printf "\n=== X_0^%o(%o)*  genus %o\n", D, N, X`g;
    ell := NumberOfEllipticPointsByCMOrder(X);
    classes := [* *]; bad := false;
    for q in Sort(SetToSequence(Keys(ell))) do
        for d in Sort(SetToSequence(Keys(ell[q]))) do
            n := ell[q][d]; if n eq 0 then continue; end if;
            ok := true; K := Rationals();
            try K := AbsoluteField(FieldOfDefinitionOfCMPoint(X, d)); catch e ok := false; end try;
            if not ok then bad := true;
               printf "   no field of definition at discriminant %o (%o points)\n", d, n; continue;
            end if;
            if Degree(K) ne n then
               printf "   WARNING d = %o: %o points but field degree %o\n", d, n, Degree(K);
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
        printf "   class: stabiliser order %o, field degree %o, %o points, discriminants %o%o\n",
           cl[1], Degree(cl[2]), cl[3], cl[4], IsOdd(cl[3]) select "   <- ODD" else "";
        if IsOdd(cl[3]) then S +:= Degree(cl[2]); Append(~odd, cl[4]); end if;
    end for;
    printf "   SUM OF FIELD DEGREES OVER ODD CLASSES = %o   (a hyperelliptic involution at "
         * "genus %o has exactly %o fixed points)\n", S, X`g, 2*X`g+2;
    printf "   VERDICT: %o%o\n",
       S gt 2*X`g+2 select "NOT HYPERELLIPTIC" else "no conclusion from this criterion",
       bad select "   [INCOMPLETE: a field was unavailable]" else "";
end for;
printf "\nELL3 DONE\n";
quit;
