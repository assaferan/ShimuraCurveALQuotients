AttachSpec("ShimuraQuotients.spec"); SetColumns(0);
curves := eval Read("data/curves_after_UpdateCurves8.dat");
printf "%o curves\n", #curves;
und := [X : X in curves | not assigned X`IsSubhyp];
sub := [X : X in curves | assigned X`IsSubhyp and X`IsSubhyp];
non := [X : X in curves | assigned X`IsSubhyp and not X`IsSubhyp];
printf "subhyperelliptic %o, not %o, undecided %o\n", #sub, #non, #und;
// genus distribution of the undecided and of the subhyperelliptic ones
gs := {X`g : X in und};
printf "undecided by genus: %o\n", [<g, #[X : X in und | X`g eq g]> : g in Sort(SetToSequence(gs))];
printf "undecided by (D,N) count: %o bases\n", #{<X`D, X`N> : X in und};
printf "undecided with N = 1: %o; D*N max %o\n", #[X : X in und | X`N eq 1], Maximum([X`D*X`N : X in und]);
printf "subhyperelliptic by genus: %o\n", [<g, #[X : X in sub | X`g eq g]> : g in Sort(SetToSequence({X`g : X in sub}))];
// how many of the subhyperelliptic genus >= 2 curves have a model stored?
quit;
