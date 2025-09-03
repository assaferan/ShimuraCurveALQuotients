declare type LogSm;
declare attributes LogSm : log_coeffs;

intrinsic LogSum() -> LogSm
{.}
    ret := New(LogSm);
    ret`log_coeffs := AssociativeArray();

    return ret;
end intrinsic;

intrinsic LogSum(oo::Infty) -> LogSm
{.}
    require oo eq Infinity() : "oo must be +Infinity";
    ret := LogSum();
    ret`log_coeffs[-1] := 1;

    return ret;
end intrinsic;

intrinsic LogSum(a::FldRatElt) -> LogSm
{.}
    ret := LogSum();
    if a eq 0 then
        ret`log_coeffs[0] := 1;
        return ret;
    end if;
    require a gt 0 : "a must be non-negative";
    fac := FactorizationOfQuotient(a);
    for p in fac do
        ret +:= LogSum(p[2],p[1]);
    end for;
    return ret;
end intrinsic;

intrinsic LogSum(a::RngIntElt) -> LogSm
{.}
    return LogSum(Rationals()!a);
end intrinsic;

intrinsic LogSum(a::FldRatElt, p::RngIntElt) -> LogSm
{The element a log p.}
    ret := New(LogSm);
    ret`log_coeffs := AssociativeArray();
    ret`log_coeffs[p] := a;

    return ret;
end intrinsic;

intrinsic LogSum(a::RngIntElt, p::RngIntElt) -> LogSm
{The element a log p.}
    return LogSum(Rationals()!a,p);
end intrinsic;

intrinsic LogSum(primes_coeffs::SetEnum[Tup]) -> LogSm
{The sum of a log p for a set of pairs <p,a>.}
    if IsEmpty(primes_coeffs) then return LogSum(); end if;
    return &+[LogSum(pa[2], pa[1]) : pa in primes_coeffs];
end intrinsic;

procedure reduce(s)
    zero_keys := [p : p in Keys(s`log_coeffs) | s`log_coeffs[p] eq 0];
    for p in zero_keys do
        Remove(~(s`log_coeffs), p);
    end for;
    return;
end procedure;

intrinsic '+'(s1::LogSm, s2::LogSm) -> LogSm
{.}
    s := New(LogSm);
    s`log_coeffs := AssociativeArray();
    for special in [0, -1] do
        if IsDefined(s1`log_coeffs, special) then
            s`log_coeffs := s1`log_coeffs;
            return s;
        end if;
        if IsDefined(s2`log_coeffs, special) then
            s`log_coeffs := s2`log_coeffs;
            return s;
        end if;
    end for;
    for p in Keys(s1`log_coeffs) do
        s`log_coeffs[p] := s1`log_coeffs[p];
    end for;
    for p in Keys(s2`log_coeffs) do
        if not IsDefined(s`log_coeffs,p) then
            s`log_coeffs[p] := 0;
        end if;
        s`log_coeffs[p] +:= s2`log_coeffs[p];
    end for;
    reduce(s);
    return s;
end intrinsic;

intrinsic '*'(a::FldRatElt, s::LogSm) -> LogSm
{.}
    s_a := New(LogSm);
    s_a`log_coeffs := AssociativeArray();
    for p in Keys(s`log_coeffs) do
        s_a`log_coeffs[p] := a*s`log_coeffs[p];
    end for;
    
    reduce(s_a);
    return s_a;
end intrinsic;

intrinsic '*'(a::RngIntElt, s::LogSm) -> LogSm
{.}
   return Rationals()!a*s;
end intrinsic;


intrinsic '-'(s1::LogSm, s2::LogSm) -> LogSm
{.}
    return s1 + (-1)*s2;
end intrinsic;

intrinsic Print(s::LogSm)
{.}
    primes := Sort([k : k in Keys(s`log_coeffs)]);
    if IsEmpty(primes) then printf "0"; return; end if;
    for j->p in primes do
        if (j ne 1) and ( s`log_coeffs[p] gt 0) then
            printf "+";
        end if;
        coeff := Sprintf("%o", s`log_coeffs[p]);
        if (Abs(s`log_coeffs[p]) eq 1) then
            coeff := coeff[1..1];
            if (j eq 1) and (s`log_coeffs[p] eq 1) then coeff := ""; end if;
        end if;
        printf "%oLog%o", coeff, p;
    end for;
    return;
end intrinsic;

intrinsic RationalNumber(s::LogSm) -> FldRatElt
{.}
    require &and[IsIntegral(coeff) : coeff in s`log_coeffs] : "s does not represent a rational number!";
    ret := &*[Rationals() | p^(Integers()!s`log_coeffs[p]) : p in Keys(s`log_coeffs)];
    if (ret eq -1) then return Infinity(); end if;
    return ret;
end intrinsic;

intrinsic IsZero(s::LogSm) -> BoolElt
{.}
    reduce(s);
    return IsEmpty(Keys(s`log_coeffs));
end intrinsic;

intrinsic 'eq'(s1::LogSm, s2::LogSm) -> BoolElt
{.}
    for special in [0,-1] do
        if special in Keys(s1`log_coeffs) or special in Keys(s2`log_coeffs) then
            return special in Keys(s2`log_coeffs) and special in Keys(s1`log_coeffs);
        end if;
    end for;
    return IsZero(s1-s2);
end intrinsic;

intrinsic SquareFree(s::LogSm) -> LogSm
{.}
    ret := LogSum();
    for p in Keys(s`log_coeffs) do
        ret`log_coeffs[p] := s`log_coeffs[p] - 2*Floor(s`log_coeffs[p]/2);
    end for;
    return ret;
end intrinsic;

intrinsic IsSquare(s::LogSm) -> BoolElt
{.}
    return SquareFree(s) eq LogSum();
end intrinsic;