// The extended dual-graph filter on the undecided curves of the final data (all with D > 1).
AttachSpec("ShimuraQuotients.spec"); SetColumns(0);
curves := eval Read("data/curves_after_UpdateCurves8.dat");
und := [i : i in [1..#curves] | not assigned curves[i]`IsSubhyp];
sub := [curves[i] : i in und];
printf "%o undecided, %o with N not squarefree\n", #sub, #[X : X in sub | not IsSquarefree(X`N)];
t0 := Realtime();
FilterByDualGraph(~sub);
printf "took %o s\n", Realtime() - t0;
dec := [X : X in sub | assigned X`IsSubhyp];
printf "decided %o: %o\n", #dec, [<X`D, X`N, Sort(SetToSequence(X`W)), X`g, X`TestInWhichProved> : X in dec];
printf "still undecided: %o (squarefree N: %o)\n", #sub - #dec, #[X : X in sub | not assigned X`IsSubhyp and IsSquarefree(X`N)];
printf "DGRUN_DONE\n";
quit;
