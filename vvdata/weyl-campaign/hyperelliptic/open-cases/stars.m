AttachSpec("ShimuraQuotients.spec"); SetColumns(0);
cs := eval Read("data/curves_after_UpdateCurves8.dat");
metric := {<10,21,{1,2,7,14}>,<10,21,{1,7,10,70}>,<14,15,{1,3,7,21}>,<14,15,{1,3,70,210}>,<14,15,{1,5,7,35}>,<14,33,{1,3,7,21,22,66,154,462}>,<15,14,{1,2,7,14}>,<22,51,{1,2,33,51,66,102,187,374}>,<22,51,{1,3,11,33,34,102,374,1122}>,<22,63,{1,2,7,9,11,14,18,22,63,77,99,126,154,198,693,1386}>,<26,85,{1,2,5,10,13,17,26,34,65,85,130,170,221,442,1105,2210}>,<34,21,{1,2,7,14,51,102,357,714}>,<34,25,{1,2,17,25,34,50,425,850}>,<38,15,{1,3,5,15,38,114,190,570}>,<57,5,{1,3,19,57}>,<62,5,{1,5,31,155}>,<65,2,{1,2,13,26}>,<86,15,{1,2,3,5,6,10,15,30,43,86,129,215,258,430,645,1290}>,<141,1,{1,3}>,<194,3,{1,6,194,291}>,<798,1,{1,3,38,114}>,<870,1,{1,5,29,145}>};
comp := {<14,51,{1,6,14,21,34,51,119,714}>,<38,15,{1,3,10,30,38,95,114,285}>,<158,3,{1,2,237,474}>};
m635 := {<6,35,{1,3,5,15}>,<6,35,{1,5,7,35}>};
rest := [X : X in cs | not assigned X`IsSubhyp and not (<X`D,X`N,X`W> in metric)
                       and not (<X`D,X`N,X`W> in comp) and not (<X`D,X`N,X`W> in m635)];
printf "remaining: %o\n", #rest;
ns := 0;
for X in rest do
    DN := X`D*X`N; full := {Integers()| d : d in Divisors(DN) | GCD(d, DN div d) eq 1};
    idx := #full div #X`W;
    if idx eq 1 then
        ns +:= 1;
        printf "STAR  X_0^%o(%o)*  g=%o  N %o\n", X`D, X`N, X`g, IsSquarefree(X`N) select "squarefree" else "NOT squarefree";
    end if;
end for;
printf "\nstar curves: %o of %o\n", ns, #rest;
idxs := {* #{Integers()| d : d in Divisors(X`D*X`N) | GCD(d, X`D*X`N div d) eq 1} div #X`W : X in rest *};
for i in Sort(SetToSequence(MultisetToSet(idxs))) do printf "index %o: %o curves\n", i, Multiplicity(idxs, i); end for;
quit;
