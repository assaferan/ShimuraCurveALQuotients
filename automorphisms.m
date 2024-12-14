function NumPointsFpd(X,p, d :cached_traces := AssociativeArray(), class_nos := AssociativeArray())
    D := X`D;
    N := X`N;
    W := X`W;
    if <D,N,2,p^d,W> notin Keys(cached_traces) then
        cached_traces[<D,N,2,p^d,W>] := TraceDNewALFixed(X`D, X`N, 2, p^d, X`W :class_nos := class_nos);
    end if;
    tpsd := cached_traces[<D,N,2,p^d,W>];

    if <D,N,2,p^(d-2),W> notin Keys(cached_traces) then 
        if d eq 1 or d eq 2 then
            cached_traces[<D,N,2,p^(d-2),W>] := 0;
        else        
            cached_traces[<D,N,2,p^(d-2),W>] := TraceDNewALFixed(X`D, X`N, 2, p^(d-2), X`W: class_nos := class_nos);
        end if;
    end if;
    tpsdmin2 := cached_traces[<D,N,2,p^(d-2),W>];

    trace_frob_n := tpsd - p*tpsdmin2;

    num_points := p^d + 1 - trace_frob_n;
    return Integers()!num_points, cached_traces;
end function;

function modr2(r)
    if r mod 2 eq 0 then 
        return 0; 
    end if;
    return 1;
end function;

function Pp(X, n, p: cached_traces := AssociativeArray(), class_nos := AssociativeArray())
    sum := 0;
    for d in Divisors(n) do
        numpts, cached_traces := NumPointsFpd(X,p,d : cached_traces := cached_traces, class_nos := class_nos);
        sum +:= modr2(MoebiusMu(n div d) * Abs(numpts));
    end for;
    return sum/n, cached_traces;
end function;

function involutioncounter(X, p,k :cached_traces := AssociativeArray(), class_nos := AssociativeArray())
    assert X`g gt 2;
    sum := 0;
    for i in [0..k] do
       valp, cached_traces:= Pp(X, 2*i+1, p :cached_traces := cached_traces, class_nos := class_nos);
     sum +:= (2*i+1)* valp;
     print "sum is", sum;
     if sum gt 2*X`g+2 then return false, sum, cached_traces; end if;
    end for;
    return true, sum, cached_traces;
end function;

