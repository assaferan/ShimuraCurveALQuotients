AttachSpec("/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/ShimuraQuotients.spec");
SetColumns(0);
load "/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/vvdata/weyl-campaign/voight-willis/vw3.m";
load "/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/vvdata/weyl-campaign/voight-willis/vwal.m";
// Each Atkin-Lehner element must act NON-TRIVIALLY on the disc and the group they generate must
// have the right order: #W* cosets, each reached.  A scalar would collapse cosets together.
for b in Split(BASES, ",") do
    xs := Split(b, "_"); D := StringToInteger(xs[1]); N := StringToInteger(xs[2]); DN := D*N;
    A := QuaternionAlgebra(D); O := QuaternionOrder(A, N);
    Wstar := [ d : d in Divisors(DN) | GCD(d, DN div d) eq 1 ];
    bad := [];
    for m in Wstar do
        if m eq 1 then continue; end if;
        ok, mu := ALElement(O, m, 3);
        if not ok then Append(~bad, <m, "not found">); continue; end if;
        scal := &and[ (A!mu)*(A!x) eq (A!x)*(A!mu) : x in Basis(O) ];
        if scal then Append(~bad, <m, "scalar">); end if;
        printf "  %o_%o m=%o: norm %o, central %o, element %o\n", D, N, m, Norm(mu), scal, mu;
    end for;
    printf "%o: %o defective of %o\n", b, #bad, #Wstar-1;
end for;
printf "ALCHK DONE\n"; quit;
