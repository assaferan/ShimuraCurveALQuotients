curves := GetHyperellipticCandidates();
unknown := [c : c in curves | not assigned c`IsSubhyp];
unknownstar := [ c : c in unknown | IsStarCurve(c)];
unknownmodstar := [ c : c in unknownstar | c`D eq 1];
lut := AssociativeArray();
for i in [1..#curves] do
    c := curves[i];
    lut[<c`D, c`N, c`W>] := i;
end for;

load "automorphisms.m";
cached_traces := AssociativeArray();
//Curves with 3 good
testcurves := [x : x in unknownstar |3 notin PrimeFactors(x`D*x`N )];
potential_curves := [];
for i->X in testcurves do
    print "doing curve number", i;
    D := X`D; N := X`N;
    p := 3; 
    b, sum := involutioncounter(X,p,3:cached_traces := cached_traces, class_nos := class_nos);
    if sum ne 0 then
        Append(~potential_curves,X);
        print "potential good one", i;
    end if;
end for;


nonhyp := [];
for i->X in potential_curves do
    print "doing curve number", i;
    D := X`D; N := X`N;
    g := X`g;
    p := 3;
    print g;
    b, sum := involutioncounter(X,p,2*g:cached_traces := cached_traces, class_nos := class_nos);
    if not b then
        print "success";
        Append(~nonhyp, X);
    end if;
end for;


//curves with 2 good
testcurves := [x : x in unknownstar |2 notin PrimeFactors(x`D*x`N )];

potential_curves := [];
for i->X in testcurves do
    print "doing curve number", i;
    D := X`D; N := X`N;
    p := 2; 
    b, sum := involutioncounter(X,p,2:cached_traces := cached_traces, class_nos := class_nos);
    if sum ne 0 then
        Append(~potential_curves,X);
        print "potential good one", i;
    end if;
end for;

nonhyp := [];
for i->X in potential_curves do
    print "doing curve number", i;
    D := X`D; N := X`N;
    g := X`g;
    print g;
    p := 2;
    b, sum := involutioncounter(X,p,2*g+2:cached_traces := cached_traces, class_nos := class_nos);
    if not b then
        print p;
        print "success";
        Append(~nonhyp, X);
    end if;
end for;
