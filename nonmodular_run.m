// run with index := {n};
// Then this will check curves[n]
AttachSpec("shimuraquots.spec");
SetVerbose("ShimuraQuotients", 0);
curves := GetHyperellipticCandidates();
// unknown := [X : X in curves | not assigned X`IsSubhyp];
X := curves[eval(index)];
is_hyp, name, fix := CheckModularNonALInvolutionTrace(X);
if is_hyp ne 0 then
    fix := -1;
end if;
if is_hyp eq -1 then
    name := "";
end if;
printf "%o:%o:%o:%o\n", index, is_hyp, name, fix;
exit;

