function NumPointsFpd(X,p, d)
    tpsd := TraceDNewALFixed(X`D, X`N, 2, p^d, X`W);
    if d eq 1 or d eq 2 then 
        tpsdmin2 := 0;
    else        
        tpsdmin2 := TraceDNewALFixed(X`D, X`N, 2, p^(d-2), X`W);
    end if;

    trace_frob_n := tpsd - p*tpsdmin2;

    num_points := p^d + 1 - trace_frob_n;
    return Integers()!num_points;
end function;

function modr2(r)
    if Integers()!r mod 2 eq 0 then 
        return 0; 
    end if;
    return 1;
end function;

function Pp(n, p)
    sum := 0;
    for d in Divisors(n) do
        sum +:= modr2(MoebiusMu(n div d) * Abs(NumPointsFpd(X,p,d)));
    end for;
    return sum/n;
end function;

function involutioncounter(X, p,k)
    assert X`g gt 2;
    print "2g+2 =", 2*X`g+2;
    sum := 0;
    for i in [0..k] do
     print "doing sum i=", i;
     sum +:= (2*i+1)* Pp(2*i+1, p);
     print sum;
     if sum gt 2*X`g+2 then return false; end if;
    end for;
    return true;
end function;

