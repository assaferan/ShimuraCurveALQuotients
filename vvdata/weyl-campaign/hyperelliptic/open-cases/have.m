// For each open genus-3 X with an elliptic quotient E = X/<W,w>: is there a stored model of E,
// and what is the genus of the star curve (a genus-0 star gives a Hauptmodul coordinate)?
AttachSpec("ShimuraQuotients.spec"); SetColumns(0);
cs := eval Read("data/curves_after_UpdateCurves8.dat");
gstar := AssociativeArray();
for X in cs do
    DN := X`D*X`N; full := {Integers()| d : d in Divisors(DN) | GCD(d, DN div d) eq 1};
    if X`W eq full then gstar[<X`D, X`N>] := X`g; end if;
end for;
open := [X : X in cs | not assigned X`IsSubhyp and X`g eq 3];
nhave := 0; ntot := 0;
for X in open do
    D := X`D; N := X`N; W := X`W; DN := D*N;
    full := {Integers()| d : d in Divisors(DN) | GCD(d, DN div d) eq 1};
    fn := Sprintf("%o/models_%o_%o.m", MD, D, N);
    hasfile := FileExists(fn);
    models := hasfile select eval (Read(fn) cat "\nreturn models;") else AssociativeArray();
    seen := {}; out := [];
    for w in full diff W do
        coset := {AtkinLehnerMul(w, u, DN) : u in W};
        if coset in seen then continue; end if; Include(~seen, coset);
        Wb := W join coset;
        if GenusShimuraCurveQuotient(D, N, Wb) ne 1 then continue; end if;
        ntot +:= 1;
        k := Sort(SetToSequence(Wb));
        got := IsDefined(models, k) and #models[k] gt 0;
        if got then nhave +:= 1; end if;
        Append(~out, Sprintf("w=%o %o", w, got select "MODEL" else "none"));
    end for;
    printf "X_0(%o,%o)/%o: file %o, g* = %o | %o\n", D, N, Sort(Setseq(W)), hasfile select "yes" else "NO",
        IsDefined(gstar, <D, N>) select gstar[<D, N>] else -1, Join(out, ", ");
end for;
printf "ROUTE2: %o of %o elliptic quotients have a stored model\n", nhave, ntot;
quit;
