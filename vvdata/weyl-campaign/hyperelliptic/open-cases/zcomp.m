// For each open curve X = X_0^D(N)/W and each odd p || N, p coprime to D: the component
// Z = X_0^D(N/p)/W' (W' = elements of W prime to p) of the special fibre at p, and its status.
AttachSpec("ShimuraQuotients.spec"); SetColumns(0);
cs := eval Read("data/curves_after_UpdateCurves8.dat");
idx := AssociativeArray();
for i->X in cs do idx[<X`D, X`N, X`W>] := i; end for;
metric := {<10,21,{1,2,7,14}>,<10,21,{1,7,10,70}>,<14,15,{1,3,7,21}>,<14,15,{1,3,70,210}>,<14,15,{1,5,7,35}>,<14,33,{1,3,7,21,22,66,154,462}>,<15,14,{1,2,7,14}>,<22,51,{1,2,33,51,66,102,187,374}>,<22,51,{1,3,11,33,34,102,374,1122}>,<22,63,{1,2,7,9,11,14,18,22,63,77,99,126,154,198,693,1386}>,<26,85,{1,2,5,10,13,17,26,34,65,85,130,170,221,442,1105,2210}>,<34,21,{1,2,7,14,51,102,357,714}>,<34,25,{1,2,17,25,34,50,425,850}>,<38,15,{1,3,5,15,38,114,190,570}>,<57,5,{1,3,19,57}>,<62,5,{1,5,31,155}>,<65,2,{1,2,13,26}>,<86,15,{1,2,3,5,6,10,15,30,43,86,129,215,258,430,645,1290}>,<141,1,{1,3}>,<194,3,{1,6,194,291}>,<798,1,{1,3,38,114}>,<870,1,{1,5,29,145}>};
open := [X : X in cs | not assigned X`IsSubhyp and not (<X`D, X`N, X`W> in metric)];
printf "%o open after #116\n", #open;
nr := 0;
for X in open do
    D := X`D; N := X`N; W := X`W;
    full := {Integers()| d : d in Divisors(D*N) | GCD(d, D*N div d) eq 1};
    line := Sprintf("X_0(%o,%o)/%o g=%o [W_full:W]=%o:", D, N, Sort(Setseq(W)), X`g, #full div #W);
    verdict := "";
    for p in PrimeDivisors(N) do
        if p eq 2 or N mod p^2 eq 0 or D mod p eq 0 then continue; end if;
        Wp := {w : w in W | w mod p ne 0};
        key := <D, N div p, Wp>;
        if not IsDefined(idx, key) then line cat:= Sprintf(" p=%o: Z not in data;", p); continue; end if;
        Z := cs[idx[key]];
        st := assigned Z`IsSubhyp select (Z`IsSubhyp select "subhyp" else "NOT hyp") else "open";
        line cat:= Sprintf(" p=%o: Z=X_0(%o,%o)/%o g=%o %o;", p, D, N div p, Sort(Setseq(Wp)), Z`g, st);
        if assigned Z`IsSubhyp and not Z`IsSubhyp and Z`g ge 3 then verdict := Sprintf(" => NOT hyperelliptic (Z at p=%o)", p); end if;
    end for;
    if verdict ne "" then nr +:= 1; end if;
    print line cat verdict;
end for;
printf "ruled out by a non-hyperelliptic component: %o of %o\n", nr, #open;
quit;
