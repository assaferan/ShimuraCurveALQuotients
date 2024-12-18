import "Caching.m" : cached_traces, SetCache, GetCache, class_nos, point_counts;



function ComputePointsViaTrace(X, p, d)
    //Compute points via trace formula on X(Fp^d) when d is between 1 and g
    g := X`g;
    D := X`D;
    N := X`N;
    W := X`W;
    assert d le g;
    b, tpsd := GetCache(<D,N,2,p^d,W>, cached_traces);
    if not b then
        tpsd := TraceDNewALFixed(D, N, 2, p^d, W );
        SetCache(<D,N,2,p^d,W>,tpsd, cached_traces);
    end if;

    if d eq 1 then
            tpsdmin2 := 0;
    else 
        b, tpsdmin2 := GetCache(<D,N,2,p^(d-2),W>, cached_traces);
        if not b then 
            tpsdmin2 := TraceDNewALFixed(D, N, 2, p^(d-2), W);  
            SetCache(<D,N,2,p^(d-2),W>,tpsdmin2, cached_traces);
        end if;
    end if;

    trace_frob_n := tpsd - p*tpsdmin2;

    num_points := p^d + 1 - trace_frob_n;
    num_points := Integers()!num_points;
    return num_points;
end function;


intrinsic NumPointsFpd(X::ShimuraQuot,p::RngIntElt, d::RngIntElt) ->RngIntElt
    {Count the number of points X(F_p^d)}
    D := X`D;
    N := X`N;
    W := X`W;
    g := X`g;

    b, Npd := GetCache(<X, p^d>, point_counts);
    if b then
        return Npd;        
    else
        Nps :=[];
        for i in [1 .. g] do
            b, Npd := GetCache(<X, p^i>, point_counts);
            if not b then
                Npd := ComputePointsViaTrace(X, p, i);
                SetCache(<X, p^i>, Npd, point_counts);
            end if;
            Append(~Nps, Npd);
        end for;
        if d in [1..g] then
            return Nps[d];
        end if;
        //now have computed all number of points up to g
        //https://arxiv.org/pdf/1704.04661
        S := [p^j +1 - Nps[j] : j in [1..g]]; 
        //S[j] := sum of alpha_i^j where alpha1, ..., alpha2g are roots of reverse P(T)
        c :=[-S[1]];
        for j in [2..g] do
            Append(~c,(-S[j] + &+[c[i]*S[j-i] : i in [1..j-1]])/j);
        end for;
        //c[i] are signed elementary symmetric polynomials of the alpha
        for j in [g+1..2*g-1] do
            Append(~c,p^(j-g)*c[2*g-j]);
        end for;
        Append(~c, p^g);

        for j in [g+1..2*g] do
            Append(~S, -(&+[c[i]*S[j-i] : i in [1..(j-1)]]+j*c[j]));
        end for;

        for j in [g+1..2*g] do
            Append(~Nps, p^j + 1 -S[j]);
            SetCache(<X, p^j>,  p^j + 1 -S[j], point_counts);
        end for; // could combine with previous loop

        if d in [g+1..2*g] then //we have computed our point count 
            return Nps[d]; 
        else //continue computing
            for j in [2*g+1, d] do
                Append(~S, -(&+[c[i]*S[j-i]: i in [1..2*g]]));
            end for;

            for j in [2*g+1, d] do
                Append(~Nps, p^j +1 -S[j]);
                SetCache(<X, p^j>,  p^j + 1 -S[j], point_counts);
            end for; //can make all this faster by not looping twice

            return Nps[d];
        end if;

    end if;

end intrinsic;

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
        sum +:= MoebiusMu(n div d) * Abs(numpts);
    end for;
    assert sum mod n eq 0;
    quot := sum div n;
    return modr2(quot);
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




intrinsic FilterStarCurvesByFpAutomorphisms(starcurves ::SeqEnum, ~curves::SeqEnum, p::RngIntElt, k::RngIntElt)
    {Choose p, if the curve X does not have Fp automorphisms, then update status in curves}

    goodredn := [x : x in starcurves |p notin PrimeFactors(x`D*x`N )];

    for i->X in goodredn do
        print "starting curve", i;
        g := X`g;
        b, sum := InvolutionCounter(X,p, k);
        if not b then
            id := X`CurveID;
            curves[id]`IsSubhyp := false;
            curves[id]`IsHyp := false;
        end if;
    end for;

end intrinsic;

