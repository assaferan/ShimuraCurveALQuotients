// Where does FilterByWeilPolynomial spend its time on a top-end curve?  Replays exactly what
// IsHypWeilPolynomial does per curve -- WeilPolynomial(X, p) for every good prime p up to the
// stage's own bound (class-number tables, budget, per-genus ceiling) -- and logs each prime's
// wall time.  Progress goes through Write() so it is not lost to stdout buffering.
// Usage: magma -b id:=13029 tag:=main weil_timing.m
AttachSpec("ShimuraQuotients.spec");
id := StringToInteger(id);
curves := eval Read("data/curves_after_UpdateCurves7.dat");
X := curves[id];
assert X`CurveID eq id;
ceiling := AssociativeArray();
ceiling[3] := 53; ceiling[4] := 53; ceiling[5] := 37; ceiling[6] := 29; ceiling[7] := 23; ceiling[8] := 17;
b  := WeilClassNumberPrimeBound(Maximum(X`W), X`g);
bb := WeilBudgetPrimeBound(Maximum(X`W), X`g);
bound := Minimum([b, bb, IsDefined(ceiling, X`g) select ceiling[X`g] else b]);
log := Sprintf("%o/weil_%o_%o.log", scratch, tag, id);
Write(log, Sprintf("curve %o: D=%o N=%o g=%o #W=%o Qmax=%o | table bound %o, budget bound %o, ceiling -> effective bound %o",
    id, X`D, X`N, X`g, #X`W, Maximum(X`W), b, bb, bound) : Overwrite := true);
t_all := Cputime(); r_all := Realtime();
for p in PrimesUpTo(bound) do
    if (X`D*X`N) mod p eq 0 then continue; end if;
    t0 := Cputime(); r0 := Realtime();
    wp := WeilPolynomial(X, p);
    Write(log, Sprintf("p=%o  cpu %o s  wall %o s  cumulative cpu %o s   wp=%o", p,
        RealField(4)!Cputime(t0), RealField(4)!(Realtime()-r0), RealField(5)!Cputime(t_all), wp));
end for;
Write(log, Sprintf("DONE  total cpu %o s  wall %o s", RealField(5)!Cputime(t_all), RealField(5)!(Realtime()-r_all)));
exit;
