potential_curves := [];
for i->X in unknownstar do
    print "doing curve number", i;
    D := X`D; N := X`N;
    primes := PrimeFactors(D*N);
    chooseprime := Set(PrimesUpTo(30)) diff Set(primes);
    p :=Sort(SetToSequence(chooseprime))[1];
    if p gt 3 then
        continue;
    else 
    b, sum := involutioncounter(X,p,5:cached_traces := cached_traces);
    if sum ne 0 then
        Append(~potential_curves,X);
        print "potential good one", i;
    end if;
    end if;
end for;

nonhyp := [];
for i->X in potential_curves do
    print "doing curve number", i;
    D := X`D; N := X`N;
    g := X`g;
    primes := PrimeFactors(D*N);
    chooseprime := Set(PrimesUpTo(30)) diff Set(primes);
    p :=Sort(SetToSequence(chooseprime))[1];
    b, sum := involutioncounter(X,p,2*g+2:cached_traces := cached_traces);
    if not b then
        print "success";
        Append(~nonhyp, X);
    end if;
end for;