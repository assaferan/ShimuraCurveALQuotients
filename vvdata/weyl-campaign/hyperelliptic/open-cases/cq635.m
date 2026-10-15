// The two open quotients of X_0^6(35) against the level-35 models of the composite-level branch
// (which rest on the composite-level m = 0 multiplier of #74).
AttachSpec("ShimuraQuotients.spec"); SetColumns(0);
function rec(D, N, W) X := CreateShimuraQuot(D, N, W); X`g := GenusShimuraCurveQuotient(D, N, W); return X; end function;
cs := [rec(6, 35, {Integers()|1,3,5,15}), rec(6, 35, {Integers()|1,5,7,35})];
FilterByCanonicalQuadrics(~cs : ModelDir := MD);
for X in cs do printf "X_0(6,35)/%o g=%o: %o %o\n", Sort(Setseq(X`W)), X`g, assigned X`IsSubhyp select (X`IsSubhyp select "HYPERELLIPTIC" else "NOT hyperelliptic") else "undecided", assigned X`TestInWhichProved select X`TestInWhichProved else ""; end for;
quit;
