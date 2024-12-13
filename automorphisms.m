// cached_traces := AssociativeArray();

function NumPointsFpd(X,p, d :cached_traces := AssociativeArray())
    D := X`D;
    N := X`N;
    W := X`W;
    if <D,N,2,p^d,W> notin Keys(cached_traces) then
        cached_traces[<D,N,2,p^d,W>] := TraceDNewALFixed(X`D, X`N, 2, p^d, X`W);
    end if;
    tpsd := cached_traces[<D,N,2,p^d,W>];

    if <D,N,2,p^(d-2),W> notin Keys(cached_traces) then 
        if d eq 1 or d eq 2 then
            cached_traces[<D,N,2,p^(d-2),W>] := 0;
        else        
            cached_traces[<D,N,2,p^(d-2),W>] := TraceDNewALFixed(X`D, X`N, 2, p^(d-2), X`W);
        end if;
    end if;
    tpsdmin2 := cached_traces[<D,N,2,p^(d-2),W>];

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

function Pp(X, n, p: cached_traces := AssociativeArray())
    sum := 0;
    for d in Divisors(n) do
        sum +:= modr2(MoebiusMu(n div d) * Abs(NumPointsFpd(X,p,d : cached_traces := cached_traces)));
    end for;
    return sum/n;
end function;

function involutioncounter(X, p,k :cached_traces := AssociativeArray())
    assert X`g gt 2;
    sum := 0;
    for i in [0..k] do
     sum +:= (2*i+1)* Pp(X, 2*i+1, p :cached_traces := cached_traces);
     if sum gt 2*X`g+2 then return false, sum; end if;
    end for;
    return true, sum;
end function;

