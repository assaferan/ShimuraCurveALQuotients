import "Caching.m" : cached_traces, SetCache, GetCache, class_nos, point_counts;

intrinsic ComputePointsViaTrace(X::ShimuraQuot, p::RngIntElt, d::RngIntElt) -> RngIntElt
    {Compute points via trace formula on X(Fp^d) (this is better if d < g) }
    g := X`g;
    D := X`D;
    N := X`N;
    W := X`W;
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
end intrinsic;


intrinsic NumPointsFpd(X::ShimuraQuot,p::RngIntElt, d::RngIntElt) ->RngIntElt
    {Count the number of points X(F_p^d) this extrapolates the trace for d >g }
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

intrinsic FilterStarCurvesByFpAutomorphisms(~curves ::SeqEnum : B:=10, k:=20)
    {Choose p leq B prime bound, if the curve X does not have Fp automorphisms, then update status in curves. Input should be starcurves.}
    ps := PrimesUpTo(B);
    for p in ps do
        for i->X in curves do
	    if assigned X`IsSubhyp then continue; end if;
            // Only works for good reduction primes
            if (X`D*X`N) mod p eq 0 then continue; end if;
            vprint ShimuraQuotients, 2: "starting curve", i;
            g := X`g;
            b, sum := InvolutionCounter(X,p, k);
            if not b then
                curves[i]`IsSubhyp := false;
                curves[i]`IsHyp := false;
                curves[i]`TestInWhichProved := Sprintf("FpAutomorphisms with p = %o, k = %o", p, k);
            end if;
        end for;
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

    for i->X in curves do
        if p in PrimeFactors(X`D*X`N) then continue; end if;
        vprint ShimuraQuotients, 2: "starting curve", i;
        b := CheckWeilPolyg3(X,p);
        if not b then
            curves[i]`IsSubhyp := false;
            curves[i]`IsHyp := false;
            curves[i]`TestInWhichProved := Sprintf("WeilPolynomialg3 with p = %o", p);
        end if;
    end for;

end intrinsic;

intrinsic PointCountParity(X::ShimuraQuot, p::RngIntElt) ->BoolElt
    {Check if #Ramification points defined over C(F_p^d) for d odd, d < 2g +2 is odd}
    g :=X`g;
    sum := 0;
    for d in [1..2*g+1] do
        if IsEven(d) then continue; end if;
        if d eq 1 then
            pts := NumPointsFpd(X, p, d);
        else
            r := Maximum(Remove(Divisors(d), Index(Divisors(d),d)));
            pts := NumPointsFpd(X, p, d) - NumPointsFpd(X,p,r);
        end if;
        // print "d is", d;
        // print GF(2)!pts;
        // if IsEven(d) then assert GF(2)!pts eq 0; end if;
        sum +:=GF(2)!(pts); 
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




intrinsic WeilClassNumberPrimeBound(Qmax::RngIntElt, g::RngIntElt) -> RngIntElt
{The largest integer b such that the Weil-polynomial trace formula at a prime p on a
genus-g Atkin-Lehner quotient whose largest involution is Qmax only requests Hurwitz
class numbers of |discriminant| <= 4*Qmax*p^g that lie within the precomputed
class-group tables (|d| < 2^40), for all p <= b. Primes p <= b are therefore served
entirely from the tables. Returns 0 if even p = 2 would exceed the tables.}
    DB := ClassNumberDataMaxAbsDisc();
    cap := DB div (4*Qmax);              // need b^g <= cap
    if cap lt 1 then return 0; end if;
    b := Iroot(cap, g);
    while b ge 1 and 4*Qmax*b^g gt DB do b -:= 1; end while;   // correct any Iroot rounding
    while 4*Qmax*(b+1)^g le DB do b +:= 1; end while;          // make maximal
    return b;
end intrinsic;

// Largest integer prime bound used for curve c: the database-tight bound, optionally
// lowered by a per-genus practical ceiling (to keep runtime reasonable at low genus,
// where the tables would otherwise permit hundreds of primes).
function EffectiveWeilPrimeBound(c, ceiling)
    b := WeilClassNumberPrimeBound(Maximum(c`W), c`g);
    if IsDefined(ceiling, c`g) and ceiling[c`g] lt b then
        b := ceiling[c`g];
    end if;
    return b;
end function;

intrinsic IsHypWeilPolynomial(X::ShimuraQuot, possible_wps ::Assoc, poss_wps_at2 ::Assoc, bound::RngIntElt) -> BoolElt, RngIntElt
    {Return false if not a hyperelliptic Weil Poly at some good prime p <= bound.}
    g := X`g;
    assert g in Keys(possible_wps);
    if g notin [3,4,5,6] then //first check at 2 by f-rank
        if (2 le bound) and (2 notin PrimeDivisors(X`D*X`N)) then
            wp := WeilPolynomial(X,2);
            slopes := SlopesWithMultiplicities(NewtonPolygon(wp,2));
            f := [i[2] : i in slopes | i[1] eq 0][1]; //find multiplicity of 0
            u := Universe(poss_wps_at2[f]);
            wp := Reverse(Coefficients(wp));
            if u!wp notin poss_wps_at2[f] then
                vprint ShimuraQuotients, 2 : u!wp;
                return false, 2;
            end if;
        end if;
    end if;

    primes := Sort([p : p in Keys(possible_wps[g]) | p le bound]);
    for p in primes do
        if p in PrimeDivisors(X`D*X`N) then continue; end if;
        wp := Reverse(Coefficients(WeilPolynomial(X,p)));
        u:= Universe(possible_wps[g][p]);
        if u!wp notin possible_wps[g][p] then
            vprint ShimuraQuotients, 2 : u!wp;
            return false, p;
        end if;
    end for;
    return true, _;

end intrinsic;


function LMFDBweilpolys(g,p)
    f :=  Read(Sprintf("data/hypg%oq%o.txt",g,p));
    lines := Split(f, "\n");
    possible_polys := { eval c : c in lines};
    return possible_polys;
end function;

function createpossiblepolys(genera, genus_bounds)
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
        away := Set(HyperellipticWeilPolysAwayFromTwo(g)); // independent of p, compute once
        for p in PrimesUpTo(genus_bounds[g]) do
            if (g eq 3 and p lt 25) or (g eq 4 and p le 5) or (g in [5,6] and p eq 2) then
                possible_wps[g][p] := LMFDBweilpolys(g,p);
            elif p ne 2 then
                possible_wps[g][p] := away;
            end if;
        end for;
    end for;
    return possible_wps, poss_wps_at2;

end function;


intrinsic FilterByWeilPolynomial(~curves::SeqEnum : genera := { c`g : c in curves | not assigned c`IsSubhyp }, prime_ceiling := AssociativeArray())
    {Filter by constraints on Weil polynomials coming from LMFDB. Each curve is checked
    at every good prime up to its own bound: the largest prime keeping all trace-formula
    class numbers inside the precomputed class-group tables (|d| < 2^40), optionally
    further lowered by the per-genus practical ceiling prime_ceiling[g]. With no
    prime_ceiling supplied a default ceiling is used (g=3,4->53, 5->37, 6->29, 7->23,
    8->17; higher genera are bounded by the tables alone).}
    if IsEmpty(genera) then return; end if;
    if IsEmpty(Keys(prime_ceiling)) then
        prime_ceiling[3] := 53; prime_ceiling[4] := 53; prime_ceiling[5] := 37;
        prime_ceiling[6] := 29; prime_ceiling[7] := 23; prime_ceiling[8] := 17;
    end if;
    // Precompute possible Weil polynomials up to the largest per-curve bound per genus.
    precompute_bd := AssociativeArray();
    for g in genera do precompute_bd[g] := 0; end for;
    for c in curves do
        if assigned c`IsSubhyp or c`g notin genera then continue; end if;
        eb := EffectiveWeilPrimeBound(c, prime_ceiling);
        if eb gt precompute_bd[c`g] then precompute_bd[c`g] := eb; end if;
    end for;
    possible_wps, poss_wps_at2 := createpossiblepolys(genera, precompute_bd);
    for i->c in curves do
        if i mod 10 eq 0 then
            vprint ShimuraQuotients, 2: i;
        end if;
        if assigned c`IsSubhyp then continue; end if;
        if c`g notin genera then continue; end if;
        b, p := IsHypWeilPolynomial(c, possible_wps, poss_wps_at2, EffectiveWeilPrimeBound(c, prime_ceiling));
        if not b then
            curves[i]`IsSubhyp := false;
            curves[i]`IsHyp := false;
            curves[i]`TestInWhichProved := Sprintf("WeilPolynomial with p = %o", p);
        end if;
    end for;
end intrinsic;


intrinsic FilterByWeilPolynomialGenusScaled(~curves::SeqEnum)
    {Filter by Weil polynomials with tight per-curve prime bounds. Each curve uses every
    good prime up to the largest value that keeps all trace-formula class numbers within
    the precomputed class-group tables (|d| < 2^40), capped by per-genus practical
    ceilings (g=3,4->53, 5->37, 6->29, 7->23, 8->17; higher genera tables-bounded).}
    FilterByWeilPolynomial(~curves);
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
{Check whether the degeneracy morphism X0(D,N) -> X0(D,N/4) induces a degree 3 map
    X0(D,N)/W -> X0(D,N/4)/W to a genus 0 curve, proving the curve is trigonal (and therefore 
    not hyperelliptic).}

    for i->c in curves do
        if assigned c`IsSubhyp then continue; end if;
        if c`N mod 4 eq 0 and IsOdd(c`N div 4) and 4 in c`W and c`g gt 2 then
            newW := {w : w in c`W | w mod 4 ne 0};
            b := CheckTrigonalByDegeneracy(c);
            if not b then //trigonal
                curves[i]`IsSubhyp := false;
                curves[i]`IsHyp := false;
                curves[i]`TestInWhichProved := Sprintf("DegeneracyMorphism to X0(%o, %o)/W is trigonal", c`D, c`N div 4);
            end if;
        end if;
    end for;
end intrinsic;

