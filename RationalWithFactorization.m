declare type FldRatFac;

declare attributes FldRatFac: Value, Fac, Sign, IsInfinity;

intrinsic Parent(x::FldRatFac) -> PowStr
{.}
    return PowerStructure(FldRatFac);
end intrinsic;

intrinsic Factorization(x::FldRatElt) -> SeqEnum
{.}
    fac_num := Factorization(Numerator(x));
    fac_denom := Factorization(Denominator(x));
    fac := fac_num cat [<fa[1], -fa[2]> : fa in fac_denom];
    return Sort(fac);
end intrinsic;

intrinsic RationalWithFactorization(x::RngIntElt) -> FldRatFac
{.}
    return RationalWithFactorization(Rationals()!x);
end intrinsic;

intrinsic RationalWithFactorization(x::FldRatElt) -> FldRatFac
{.}
    ret := New(FldRatFac);
    ret`Value := x;
    if (x eq 0) then 
        ret`Fac := [];
    else
        ret`Fac := Factorization(x);
    end if;
    ret`Sign := Sign(x);
    ret`IsInfinity := false;
    return ret;
end intrinsic;

intrinsic RationalWithFactorization(x::Infty) -> FldRatFac
{.}
    ret := New(FldRatFac);
    ret`Value := x;
    ret`Fac := [];
    ret`Sign := Sign(x);
    ret`IsInfinity := true;
    return ret;
end intrinsic;


intrinsic Print(x::FldRatFac, level::MonStgElt)
{.}
    if level eq "Magma" then 
        printf "RationalWithFactorization(%o)", x`Value;
        return;
    end if;
    if x`Sign eq 0 then printf "0"; return; end if;
    if x`Sign eq -1 then printf "-"; end if;
    if x`IsInfinity then printf "oo"; return; end if;
    if IsEmpty(x`Fac) then printf "1"; return; end if;
    printf "%o^%o", x`Fac[1][1], x`Fac[1][2];
    for fa in x`Fac[2..#x`Fac] do
        printf " * %o^%o", fa[1], fa[2];
    end for;
    return;
end intrinsic;