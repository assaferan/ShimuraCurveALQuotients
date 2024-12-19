import "Caching.m" : cached_traces, SetCache, GetCache, class_nos, point_counts;

function ComputePointsViaTrace(X, p, d)
    //Compute points via trace formula on X(Fp^d) when d is between 1 and g
    g := X`g;
    D := X`D;
    N := X`N;
    W := X`W;
    // assert d le g;
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

    assert D mod p ne 0;
    assert N mod p ne 0;
    assert g gt 0;

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
            Append(~c,-(S[j] + &+[c[i]*S[j-i] : i in [1..j-1]]) /j);
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
            for j in [2*g+1.. d] do
                Append(~S, -(&+[c[i]*S[j-i]: i in [1..2*g]]));
            end for;

            for j in [2*g+1..d] do
                Append(~Nps, p^j +1 -S[j]);
                SetCache(<X, p^j>,  p^j + 1 -S[j], point_counts);
            end for; //can make all this faster by not looping twice

            return Nps[d];
        end if;

    end if;

end intrinsic;

intrinsic WeilPolynomial(X::ShimuraQuot,p::RngIntElt) -> RngUPolElt
    {Return Weil polynomial}
    D := X`D;
    N := X`N;
    W := X`W;
    g := X`g;

    assert D mod p ne 0;
    assert N mod p ne 0;
    assert g gt 0;
    Nps :=[];
    for i in [1 .. g] do
        b, Npd := GetCache(<X, p^i>, point_counts);
        if not b then
            Npd := ComputePointsViaTrace(X, p, i);
            SetCache(<X, p^i>, Npd, point_counts);
        end if;
        Append(~Nps, Npd);
    end for;

    //now have computed all number of points up to g
    S := [p^j +1 - Nps[j] : j in [1..g]]; 
    //S[j] := sum of alpha_i^j where alpha1, ..., alpha2g are roots of reverse P(T)
    c :=[-S[1]];
    for j in [2..g] do
        Append(~c,-(S[j] + &+[c[i]*S[j-i] : i in [1..j-1]]) /j);
    end for;
    //c[i] are signed elementary symmetric polynomials of the alpha
    for j in [g+1..2*g-1] do
        Append(~c,p^(j-g)*c[2*g-j]);
    end for;
    Append(~c, p^g);

    _<T>:=PolynomialRing(Integers());
    wp := &+[c[i]*T^(2*g-i) : i in [1..2*g]] + T^(2*g);
    return wp;

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
     vprint ShimuraQuotients, 3: "sum is", sum;
     if sum gt 2*X`g+2 then return false, sum; end if;
    end for;
    return true, sum;
end intrinsic;

intrinsic FilterStarCurvesByFpAutomorphisms(starcurves ::SeqEnum, ~curves::SeqEnum, p::RngIntElt, k::RngIntElt)
    {Choose p, if the curve X does not have Fp automorphisms, then update status in curves}

    goodredn := [x : x in starcurves |p notin PrimeFactors(x`D*x`N )];

    for i->X in goodredn do
        vprint ShimuraQuotients, 2: "starting curve", i;
        g := X`g;
        b, sum := InvolutionCounter(X,p, k);
        if not b then
            id := X`CurveID;
            curves[id]`IsSubhyp := false;
            curves[id]`IsHyp := false;
        end if;
    end for;
end intrinsic;

intrinsic CheckWeilPolyg3(X::ShimuraQuot, p::RngIntElt)->BoolElt
    {Apply Theorem 2.8 of arxiv.org/pdf/2002.02067 
    return false if not hyperelliptic, otherwise return true (unknown)}
    g:= X`g;
    assert g eq 3;
    assert p ne 2;
    poly := WeilPolynomial(X,p);
    if Integers()!(Coefficient(poly, 2)/p) mod 2 eq 0 and Coefficient(poly,3) mod 2 eq 1 then
        return false;
    else
        return true;
    end if;
end intrinsic;

intrinsic FilterByWeilPolynomialg3(~curves::SeqEnum, p)
    {This is made redundant in small p by later data}

    goodredn := [x : x in curves |p notin PrimeFactors(x`D*x`N )];

    for i->X in goodredn do
        vprint ShimuraQuotients, 2: "starting curve", i;
        g := X`g;
        b := CheckWeilPolyg3(X,p);
        if not b then
            id := X`CurveID;
            curves[id]`IsSubhyp := false;
            curves[id]`IsHyp := false;
        end if;
    end for;

end intrinsic;

intrinsic PointCountParity(X::ShimuraQuot, p::RngIntElt) ->BoolElt
    {Check if #Ramification points defined over C(F_p^d) for d odd, d < 2g +2 is odd}
    g :=X`g;
    sum := GF(2)!0;
    for d in [1..2*g+1] do
        if d eq 1 then
            pts := NumPointsFpd(X, p, d);
        else
            pts := NumPointsFpd(X, p, d) - NumPointsFpd(X,p,d-1);
        end if;
        sum +:=GF(2)!sum; 
    end for;
    if sum eq GF(2)!1 then
        return false; //return false if not hyperelliptic
    else
        return true; //unclear if hyperelliptic or not
    end if;
end intrinsic;

intrinsic HyperellipticWeilPolysAtTwo(f::RngIntElt) -> SeqEnum
{Returns all mod 2 classes of hyperelliptic Weil polynomials over F_2 of 2-rank f.}
    ds := RestrictedPartitions(f+1, {j : j in [1..f+1] | IsOdd(j)});
    _<t> := PolynomialRing(GF(2));
    ret := [];
    for d in ds do
        act_W := &*[t^x - 1 : x in d];
        Append(~ret, act_W div (t-1));
    end for;
    return ret;
end intrinsic;

intrinsic HyperellipticWeilPolysAwayFromTwo(g::RngIntElt) -> SeqEnum
{Returns all mod 2 classes of Weil polynomials over F_p (p odd) of a hyperelliptic curve of genus g.}
    ds := RestrictedPartitions(2*g+2, {j : j in [1..2*g+2] | IsOdd(j)});
    _<t> := PolynomialRing(GF(2));
    ret := [];
    for d in ds do
        act_W := &*[t^x - 1 : x in d];
        Append(~ret, act_W div (t-1)^2);
    end for;
    return ret;
end intrinsic;

intrinsic HypWeilPolynomialG3(X::ShimuraQuot) -> BoolElt
    {check against LMFDB data for g=3}
    assert X`g eq 3;

    for p in [2, 3, 5, 7, 11, 13, 17, 19, 23] do
        if p in PrimeDivisors(X`D*X`N) then 
            continue;
        end if;
        f :=  Read("nonhyp" cat Sprint(p) cat ".txt");
        lines := Split(f, "\n");
        possible_polys := { eval c : c in lines};

        poly := WeilPolynomial(X, p);
        wp := Reverse(Coefficients(poly));
        if wp in possible_polys then
            vprint  ShimuraQuotients, 3: wp;
            return false;
        end if;

    end for;
    return true;

end intrinsic;


intrinsic HypWeilPolynomialG4(X::ShimuraQuot) -> BoolElt
    {check against LMFDB data for g=4}
    assert X`g eq 4;

    for p in [2, 3, 5] do
        if p in PrimeDivisors(X`D*X`N) then 
            continue;
        end if;
        f :=  Read("nonhypg4q" cat Sprint(p) cat ".txt");
        lines := Split(f, "\n");
        possible_polys := { eval c : c in lines};

        poly := WeilPolynomial(X, p);
        wp := Reverse(Coefficients(poly));
        if wp in possible_polys then
            vprint  ShimuraQuotients, 3: wp;
            return false;
        end if;

    end for;
    return true;

end intrinsic;

intrinsic FilterByWeilPolynomialG3G4(~curves::SeqEnum)
    {Filter by constraints on weil polynomials coming from LMFDB}
    for i->c in curves do
        if i mod 10 eq 0 then
            vprint ShimuraQuotients, 2: i;
        end if;
        if assigned c`IsSubhyp then continue; end if;
        g:= c`g;
        if g eq 3 then
            b := HypWeilPolynomialG3(c);
            if not b then
                curves[i]`IsSubhyp := false;
                curves[i]`IsHyp := false;
            end if;
        elif g eq 4 then
            b := HypWeilPolynomialG4(c);
            if not b then
                curves[i]`IsSubhyp := false;
                curves[i]`IsHyp := false;
            end if;
        else
            continue;
        end if;
    end for;
end intrinsic;

