// Find a genus-3 Atkin-Lehner quotient of some X_0^D(1) whose hyperellipticity is already
// decided, to serve as the control for the canonical-quadric count built on Nelson's evaluator.
// N = 1 so that the D-new space of level D consists of newforms of level exactly D and nothing
// is old; Theorem 3.1 is stated for a newform of level d_B N.
//
//   magma -b DMAX:=400 scan3.m < /dev/null
AttachSpec("/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/ShimuraQuotients.spec");
SetColumns(0);
dmax := StringToInteger(DMAX);
MD := "/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/rerun/data/models";
for D in [1..dmax] do
    f := Factorisation(D);
    if #f eq 0 or not IsSquarefree(D) or IsOdd(#f) then continue; end if;
    ds := Divisors(D);
    for W in Subsets(Set(ds)) do
        if #W eq 0 or not (1 in W) then continue; end if;
        // W must be a subgroup of the Atkin-Lehner group (Hall divisors under m*n/gcd^2)
        ok := true;
        for a in W do for b in W do
            c := (a*b) div GCD(a,b)^2;
            if not c in W then ok := false; end if;
        end for; end for;
        if not ok then continue; end if;
        g := GenusShimuraCurveQuotient(D, 1, W);
        if g ne 3 then continue; end if;
        fn := MD cat "/models_" cat IntegerToString(D) cat "_1.m";
        has := System("test -f " cat fn) eq 0;
        printf "D=%-5o W=%-28o genus 3   model file %o\n", D, Sort(SetToSequence(W)),
               has select "PRESENT" else "absent";
    end for;
end for;
printf "SCAN3 DONE\n";
quit;
