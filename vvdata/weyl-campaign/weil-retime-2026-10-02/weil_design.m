// Pick curves for the Weil-stage re-timing: one per (genus, #W, level bucket) cell where one
// exists, plus the five already timed.  Prints id, shape, prime bound, admissible primes and the
// trace-formula term count T = sum_p sum_{w in W} 2 sqrt(p^g / Q_w) at n = p^g.
AttachSpec("ShimuraQuotients.spec");
curves := eval Read("data/curves_after_UpdateCurves7.dat");
ceiling := AssociativeArray();
ceiling[3] := 53; ceiling[4] := 53; ceiling[5] := 37; ceiling[6] := 29; ceiling[7] := 23; ceiling[8] := 17;
R := RealField(6);
function info(X)
    b := WeilClassNumberPrimeBound(Maximum(X`W), X`g); bb := WeilBudgetPrimeBound(Maximum(X`W), X`g);
    bound := Minimum([b, bb, IsDefined(ceiling, X`g) select ceiling[X`g] else b]);
    ps := [p : p in PrimesUpTo(bound) | (X`D*X`N) mod p ne 0];
    T := &+[R | 2 * Sqrt(R!p^X`g / R!Q) : p in ps, Q in X`W];
    return bound, ps, T;
end function;
bucket := func<DN | DN lt 500 select "S" else (DN lt 3000 select "M" else "L")>;
chosen := AssociativeArray();
fixed := [1071, 7296, 13029, 9755, 785];
for id in fixed do chosen[id] := true; end for;
for X in curves do
    if X`g lt 3 or X`g gt 7 then continue; end if;
    key := <X`g, #X`W, bucket(X`D*X`N)>;
    if exists{c : c in Keys(chosen) | <curves[c]`g, #curves[c]`W, bucket(curves[c]`D*curves[c]`N)> eq key} then continue; end if;
    // prefer a curve with many admissible primes (the stage's full work), and avoid level 1 oddities
    chosen[X`CurveID] := true;
end for;
ids := Sort([id : id in Keys(chosen)]);
printf "%-6o %-10o %-2o %-3o %-6o %-6o %-5o %-24o %o\n", "id", "D_N", "g", "#W", "Qmax", "DN", "bound", "primes", "terms";
tot := R!0;
for id in ids do
    X := curves[id]; bound, ps, T := info(X); tot +:= T;
    printf "%-6o %-10o %-2o %-3o %-6o %-6o %-5o %-24o %o\n", id, Sprintf("%o_%o", X`D, X`N), X`g, #X`W, Maximum(X`W), X`D*X`N, bound, ps, Floor(T);
end for;
printf "%o curves, total term count %o\n", #ids, Floor(tot);
