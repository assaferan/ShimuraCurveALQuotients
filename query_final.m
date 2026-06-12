SetQuitOnError(true);
AttachSpec("ShimuraQuotients.spec");
Attach("GeneralizedComplicatedFixedPoints.m");
SetVerbose("ShimuraQuotients", 0);

curves := eval Read("data/par_run/curves_after_UpdateCurves8.dat");
n := #curves;

// ---- Q3: classification census ----
subhyp := 0; nonsub := 0; open := 0;
p1 := 0; ec := 0; hyp := 0;
for X in curves do
    if assigned X`IsSubhyp then
        if X`IsSubhyp then subhyp +:= 1; else nonsub +:= 1; end if;
    else
        open +:= 1;
    end if;
    if assigned X`IsP1 and X`IsP1 then p1 +:= 1; end if;
    if assigned X`IsEC and X`IsEC then ec +:= 1; end if;
    if assigned X`IsHyp and X`IsHyp then hyp +:= 1; end if;
end for;
printf "TOTAL curves: %o\n", n;
printf "  classified : %o   (subhyperelliptic %o  |  non-subhyperelliptic %o)\n", subhyp+nonsub, subhyp, nonsub;
printf "    of subhyp: P1=%o  elliptic=%o  hyperelliptic=%o\n", p1, ec, hyp;
printf "  remaining (open): %o\n", open;

// ---- Q2: inspect X_0(152) family (D=1, N=152) ----
print "";
print "=== D=1, N=152 curves in the run ===";
N := 152; D := 1; DN := D*N;
printf "genus X_0(152) = %o ; star quotient genus = %o\n",
    GenusShimuraCurve(D,N), GenusShimuraCurveQuotient(D,N, {Integers()| d : d in Divisors(DN) | GCD(d, DN div d) eq 1});
for X in curves do
    if X`D eq 1 and X`N eq 152 then
        det := assigned X`IsSubhyp select (X`IsSubhyp select "SUBHYP" else "NON-subhyp") else "open";
        proved := assigned X`TestInWhichProved select X`TestInWhichProved else "-";
        printf "  W=%o g=%o  [%o]  via: %o\n", Sort(SetToSequence(X`W)), X`g, det, proved;
        // does the generalized filter fire on it?
        if X`g ge 3 then
            ok, wit := CheckGeneralizedComplicatedFixedPoints(X);
            printf "      generalized-Prop6 fires: %o%o\n", ok, ok select (" -> " cat wit) else "";
        end if;
        // geometric fixed-point counts of AL involutions w_Q (Q not in W) on this quotient
        for Q in [q : q in Divisors(DN) | GCD(q, DN div q) eq 1 and q notin X`W] do
            printf "      AL w_%o: nu_on_X=%o  count_on_quotient=%o  h(-4*%o)=%o\n",
                Q, NumFixedPoints(D,N,Q), CountFixedPointsOnQuotient(Q, X), Q, ClassNumber(-4*Q);
        end for;
    end if;
end for;
quit;
