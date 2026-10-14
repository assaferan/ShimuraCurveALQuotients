// The open cases after #113 (the data here) and #114 (the 83 the extended dual graph rules out).
AttachSpec("ShimuraQuotients.spec"); SetColumns(0);
curves := eval Read("data/curves_after_UpdateCurves8.dat");
dec := {};
for l in Split(Read(DEC), "\n") do p := Split(l, " "); Include(~dec, <StringToInteger(p[1]), StringToInteger(p[2]), {Integers()| StringToInteger(x) : x in Split(p[3], ",")}>); end for;
und := [X : X in curves | not assigned X`IsSubhyp and not (<X`D, X`N, X`W> in dec)];
printf "%o open cases\n", #und;
out := "";
for X in und do
    DN := X`D*X`N; full := {Integers()| d : d in Divisors(DN) | GCD(d, DN div d) eq 1};
    star := [Y : Y in curves | Y`D eq X`D and Y`N eq X`N and Y`W eq full];
    out cat:= Sprintf("%o %o %o g=%o index=%o star_genus=%o N_squarefree=%o\n", X`D, X`N, Sort(SetToSequence(X`W)), X`g, #full div #X`W, star[1]`g, IsSquarefree(X`N));
end for;
Write(OUT, out : Overwrite := true);
printf "by genus %o; index %o; squarefree N %o; star curves %o\n",
    [<g, #[X : X in und | X`g eq g]> : g in [3..7]],
    [<i, #[X : X in und | i eq (#{Integers()| d : d in Divisors(X`D*X`N) | GCD(d, X`D*X`N div d) eq 1}) div #X`W]> : i in [1, 2, 4]],
    #[X : X in und | IsSquarefree(X`N)], #[X : X in und | #X`W eq #{Integers()| d : d in Divisors(X`D*X`N) | GCD(d, X`D*X`N div d) eq 1}];
quit;
