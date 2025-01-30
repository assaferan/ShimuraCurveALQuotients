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
            curves[id]`TestInWhichProved := "FpAutomorphisms";
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
            curves[id]`TestInWhichProved := "WeilPolynomialg3";
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
        Append(~ret, Reverse(Coefficients(act_W div (t-1))));
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
        Append(~ret, Reverse(Coefficients(act_W div (t-1)^2)));
    end for;
    return ret;
end intrinsic;




intrinsic IsHypWeilPolynomial(X::ShimuraQuot, possible_wps ::Assoc, poss_wps_at2 ::Assoc) -> BoolElt
    {Return false if not a hyperellitic Weil Poly at some prime p}
    g := X`g;
    assert g in Keys(possible_wps);
    if g notin [3,4,5,6] then //first check at 2 by f-rank
        if 2 notin PrimeDivisors(X`D*X`N) then
            wp := WeilPolynomial(X,2);
            slopes := SlopesWithMultiplicities(NewtonPolygon(wp,2));
            f := [i[2] : i in slopes | i[1] eq 0][1]; //find multiplicity of 0
            u := Universe(poss_wps_at2[f]);
            if u!wp notin poss_wps_at2[f] then
                vprint ShimuraQuotients, 2 : u!wp;
                return false;
            end if;
        end if;
    end if;

    primes := Keys(possible_wps[g]);
    for p in primes do
        if p in PrimeDivisors(X`D*X`N) then continue; end if;
        wp := Reverse(Coefficients(WeilPolynomial(X,p)));
        u:= Universe(possible_wps[g][p]);
        if u!wp notin possible_wps[g][p] then
            vprint ShimuraQuotients, 2 : u!wp;
            return false;
        end if;
    end for;
    return true;

end intrinsic;


function LMFDBweilpolys(g,p)
    f :=  Read(Sprintf("hypg%oq%o.txt",g,p));
    lines := Split(f, "\n");
    possible_polys := { eval c : c in lines};
    return possible_polys;
end function;

function createpossiblepolys(genera : bd := 25)
    possible_wps := AssociativeArray();
    poss_wps_at2 := AssociativeArray();
    //first do 2
    largest_g := Maximum(genera);
    for f in [0..largest_g] do
        poss_wps_at2[f] := HyperellipticWeilPolysAtTwo(f);
    end for;


    for g in genera do
        assert g notin [0,1,2];
        possible_wps[g] := AssociativeArray();
        for p in PrimesUpTo(bd) do
            if (g eq 3 and p lt 25) or (g eq 4 and p le 5) or (g in [5,6] and p eq 2) then
                polys := LMFDBweilpolys(g,p);
                possible_wps[g][p] := polys;
            elif p ne 2 then
                polys := Set(HyperellipticWeilPolysAwayFromTwo(g));
                possible_wps[g][p] := polys;
            end if;
        end for;
    end for;
    return possible_wps, poss_wps_at2;

end function;


intrinsic FilterByWeilPolynomial(~curves::SeqEnum : bd := 25, genera := { c`g : c in curves | not assigned c`IsSubhyp })
    {Filter by constraints on weil polynomials coming from LMFDB}
    // genera := { c`g : c in curves | not assigned c`IsSubhyp };
    possible_wps, poss_wps_at2 := createpossiblepolys(genera :bd := bd);
    for i->c in curves do
        if i mod 10 eq 0 then
            vprint ShimuraQuotients, 2: i;
        end if;
        if assigned c`IsSubhyp then continue; end if;
        if c`g notin genera then continue; end if;
        b := IsHypWeilPolynomial(c, possible_wps, poss_wps_at2);
        if not b then
            curves[i]`IsSubhyp := false;
            curves[i]`IsHyp := false;
            curves[i]`TestInWhichProved := "WeilPolynomial";
        end if;
    end for;
end intrinsic;


intrinsic CheckTrigonalByDegeneracy(X::ShimuraQuot) -> BoolElt
    {}
    N := X`N;
    W := X`W;
    D := X`D;
    assert N mod 4 eq 0 and IsOdd(N div 4) and 4 in W and X`g gt 2;
    newW := {w : w in W | w mod 4 ne 0};
    g := GenusShimuraCurveQuotient(D,N div 4, newW);
    if g eq 0 then //trigonal
        return false;
    else 
        return true; //don't know
    end if;

end intrinsic;

intrinsic FilterByDegeneracyMorphism(~curves::SeqEnum)
    {Check whether the degeneracy morphism X0(N) -> X0(N/4) induces a degree 3 map
    X0(N)* -> X0(N/4)* to a genus 0 curve, proving the curve is trigonal (and therefore 
    not hyperelliptic).}

    for i->c in curves do
        if assigned c`IsSubhyp then continue; end if;
        if c`N mod 4 eq 0 and IsOdd(c`N div 4) and 4 in c`W and c`g gt 2 then
            newW := {w : w in c`W | w mod 4 ne 0};
            b := CheckTrigonalByDegeneracy(c);
            if not b then //trigonal
                curves[i]`IsSubhyp := false;
                curves[i]`IsHyp := false;
                curves[i]`TestInWhichProved := "DegeneracyMorphism";
            end if;
        end if;
    end for;
end intrinsic;

