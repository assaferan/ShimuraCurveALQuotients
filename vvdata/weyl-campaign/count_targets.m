SetQuitOnError(true); SetColumns(0);
AttachSpec("ShimuraQuotients.spec");
curves := GetHyperellipticCandidates();
stars := [X : X in curves | IsStarCurve(X)];
n_D1 := 0; n_nocov := 0; tgt := [];
for X in stars do
    if X`D le 1 then n_D1 +:= 1; continue; end if;
    if IsEmpty(X`CoveredBy) then n_nocov +:= 1; continue; end if;
    gs := [curves[i]`g : i in X`CoveredBy];
    dem := Maximum([2*g+5 : g in gs]);
    Append(~tgt, <X`D, X`N, dem>);
end for;
printf "stars %o | dropped D<=1: %o | dropped no covers: %o | GENUINE TARGETS: %o\n",
       #stars, n_D1, n_nocov, #tgt;
h := AssociativeArray();
for t in tgt do
    if not IsDefined(h, t[3]) then h[t[3]] := 0; end if;
    h[t[3]] +:= 1;
end for;
printf "demand histogram (2g+5 -> count):\n";
for k in Sort(SetToSequence(Keys(h))) do printf "  %o -> %o\n", k, h[k]; end for;
for t in tgt do printf "TGT %o_%o dem %o\n", t[1], t[2], t[3]; end for;
printf "COUNT_DONE\n";
quit;
