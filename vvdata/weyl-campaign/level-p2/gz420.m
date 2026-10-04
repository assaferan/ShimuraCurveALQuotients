SetColumns(0);
tab := [<-7,-7>,<-15,-5/3>,<-16,1>,<-28,1/9>,<-60,9>,<-84,-3>,<-100,1/5>,<-112,25>,<-120,-1/3>,<-148,37/9>,<-168,0>,<-228,-25/3>,<-232,-32>,<-280,-35/9>,<-312,-8/3>,<-372,-3/4>,<-408,-75>,<-532,-19/4>,<-708,25/48>,<-840,-16/3>];
for cand in [Rationals() | 21, 7/3] do
    bad := [];
    for t in tab do
        df := cand - t[2]; if df eq 0 then continue; end if;
        ps := Set(PrimeDivisors(Numerator(df))) join Set(PrimeDivisors(Denominator(df)));
        for p in ps do
            if p le 7 then continue; end if;
            M := 420*AbsoluteValue(t[1]);
            // Gross-Zagier: p can divide the difference only if 4p | 420|d| - x^2 for some x, i.e. 420|d| is a square mod p
            if M mod p eq 0 or KroneckerSymbol(M, p) eq 1 then continue; end if;
            Append(~bad, <t[1], df, p>);
        end for;
    end for;
    printf "s(-420) = %o: primes a Gross-Zagier bound forbids: %o\n", cand, bad;
end for;
