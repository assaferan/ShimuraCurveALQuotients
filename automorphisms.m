import "Caching.m" : cached_traces, SetCache, GetCache, class_nos;

function NumPointsFpd(X,p, d )
    D := X`D;
    N := X`N;
    W := X`W;
    b, tpsd := GetCache(<D,N,2,p^d,W>, cached_traces);
    if not b then
        tpsd := TraceDNewALFixed(X`D, X`N, 2, p^d, X`W );
        SetCache(<D,N,2,p^d,W>,tpsd, cached_traces);
    end if;
    


    if d eq 1 or d eq 2 then
            tpsdmin2 := 0;
    else 
        b, tpsdmin2 := GetCache(<D,N,2,p^(d-2),W>, cached_traces);
        if not b then 
            tpsdmin2 := TraceDNewALFixed(X`D, X`N, 2, p^(d-2), X`W);  
            SetCache(<D,N,2,p^(d-2),W>,tpsdmin2, cached_traces);
        end if;
    end if;

    trace_frob_n := tpsd - p*tpsdmin2;

    num_points := p^d + 1 - trace_frob_n;
    return Integers()!num_points;
end function;

function modr2(r)
    if r mod 2 eq 0 then 
        return 0; 
    end if;
    return 1;
end function;

function Pp(X, n, p)
    sum := 0;
    for d in Divisors(n) do
        numpts := NumPointsFpd(X,p,d);
        sum +:= modr2(MoebiusMu(n div d) * Abs(numpts));
    end for;
    return sum/n;
end function;

intrinsic InvolutionCounter(X ::ShimuraQuot, p ::RngIntElt,k ::RngIntElt) -> RngIntElt
    {Get upper bound on number of involutions over Fp by counting ramification points. If exceeds 2g+2 then return false.
    If false, and p of good reduction, then X not hyperelliptic (does not have involution over Fp). It is enough to choose k= 2g+2}
    assert X`g gt 2;
    sum := 0;
    for i in [0..k] do
       valp:= Pp(X, 2*i+1, p);
     sum +:= (2*i+1)* valp;
     print "sum is", sum;
     if sum gt 2*X`g+2 then return false, sum; end if;
    end for;
    return true, sum;
end intrinsic;

