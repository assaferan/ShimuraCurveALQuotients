// Enumerate the model backlog, sorted by the cost law #div(M), M = 4*D*N/2^v2(D).
// A "target" star base = D >= 2, squarefree N, and at least one subhyperelliptic quotient.
// Prints one line per base: BASE D N ndivM M nsubhyp DONE|MISSING
SetQuitOnError(true);
AttachSpec("ShimuraQuotients.spec");
SetVerbose("ShimuraQuotients", 0);
curves := eval Read("data/curves_after_UpdateCurves8.dat");
printf "LOADED %o\n", #curves;

subhyp := AssociativeArray();   // <D,N> -> count of subhyperelliptic quotients
for X in curves do
    if X`D lt 2 then continue; end if;
    if not IsSquarefree(X`N) then continue; end if;
    if not assigned X`IsSubhyp then continue; end if;
    if not X`IsSubhyp then continue; end if;
    k := <X`D, X`N>;
    if not IsDefined(subhyp, k) then subhyp[k] := 0; end if;
    subhyp[k] +:= 1;
end for;

ks := Sort([k : k in Keys(subhyp)], func<a,b | a[1]*a[2] - b[1]*b[2]>);
rows := [];
for k in ks do
    D := k[1]; N := k[2];
    M := 4 * (D*N div 2^Valuation(D,2));
    nd := #Divisors(M);
    done := "MISSING";
    try
        _ := Read(Sprintf("data/models/models_%o_%o.m", D, N));
        done := "DONE";
    catch e
        done := "MISSING";
    end try;
    Append(~rows, <nd, M, D, N, subhyp[k], done>);
end for;
Sort(~rows, func<a,b | a[1] ne b[1] select a[1]-b[1] else a[2]-b[2]>);
for r in rows do
    printf "BASE %o %o ndiv %o M %o nsub %o %o\n", r[3], r[4], r[1], r[2], r[5], r[6];
end for;
printf "TOTAL %o\n", #rows;
printf "DONEC %o\n", #[r : r in rows | r[6] eq "DONE"];
printf "MISSC %o\n", #[r : r in rows | r[6] eq "MISSING"];
printf "FINISHED\n";
quit;
