// A genus-3 curve with two distinct elliptic Atkin-Lehner quotients is the fibre product of them
// over the star curve (bielliptic construction, no Borcherds run for X itself: only the two
// genus-1 double covers y^2 = f_i(t) of the genus-0 star are needed).
// Which open curves qualify, and how far is the level from the #div(M) >= 24 Normaliz wall?
AttachSpec("ShimuraQuotients.spec"); SetColumns(0);
cs := eval Read("data/curves_after_UpdateCurves8.dat");
gstar := AssociativeArray();
for X in cs do
    DN := X`D*X`N; full := {Integers()| d : d in Divisors(DN) | GCD(d, DN div d) eq 1};
    if X`W eq full then gstar[<X`D, X`N>] := X`g; end if;
end for;
metric := {<10,21,{1,2,7,14}>,<10,21,{1,7,10,70}>,<14,15,{1,3,7,21}>,<14,15,{1,3,70,210}>,<14,15,{1,5,7,35}>,<14,33,{1,3,7,21,22,66,154,462}>,<15,14,{1,2,7,14}>,<22,51,{1,2,33,51,66,102,187,374}>,<22,51,{1,3,11,33,34,102,374,1122}>,<22,63,{1,2,7,9,11,14,18,22,63,77,99,126,154,198,693,1386}>,<26,85,{1,2,5,10,13,17,26,34,65,85,130,170,221,442,1105,2210}>,<34,21,{1,2,7,14,51,102,357,714}>,<34,25,{1,2,17,25,34,50,425,850}>,<38,15,{1,3,5,15,38,114,190,570}>,<57,5,{1,3,19,57}>,<62,5,{1,5,31,155}>,<65,2,{1,2,13,26}>,<86,15,{1,2,3,5,6,10,15,30,43,86,129,215,258,430,645,1290}>,<141,1,{1,3}>,<194,3,{1,6,194,291}>,<798,1,{1,3,38,114}>,<870,1,{1,5,29,145}>};
comp := {<14,51,{1,6,14,21,34,51,119,714}>,<38,15,{1,3,10,30,38,95,114,285}>,<158,3,{1,2,237,474}>};
m635 := {<6,35,{1,3,5,15}>,<6,35,{1,5,7,35}>};
for X in cs do
    if assigned X`IsSubhyp or X`g ne 3 then continue; end if;
    D := X`D; N := X`N; W := X`W; DN := D*N; key := <D, N, W>;
    full := {Integers()| d : d in Divisors(DN) | GCD(d, DN div d) eq 1};
    seen := {}; ell := [];
    for w in full diff W do
        coset := {AtkinLehnerMul(w, u, DN) : u in W};
        if coset in seen then continue; end if; Include(~seen, coset);
        if GenusShimuraCurveQuotient(D, N, W join coset) eq 1 then Append(~ell, w); end if;
    end for;
    if #ell lt 2 then continue; end if;
    st := key in metric select "#116" else (key in comp select "component" else (key in m635 select "6_35 models" else "OPEN"));
    printf "X_0(%o,%o)/%o: %o elliptic quotients (w = %o), g* = %o, #div(%o) = %o %o -- %o\n",
        D, N, Sort(Setseq(W)), #ell, ell, gstar[<D, N>], DN, #Divisors(DN),
        #Divisors(DN) ge 24 select "(ABOVE the wall)" else "(below the wall)", st;
end for;
quit;
