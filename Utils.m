// General utility intrinsics.

intrinsic WriteStderr(s::MonStgElt)
{ write to stderr }
  E := Open("/dev/stderr", "a");
  Write(E, s);
  Flush(E);
end intrinsic;

intrinsic WriteStderr(e::Err)
{ write to stderr }
  WriteStderr(Sprint(e) cat "\n");
end intrinsic;

intrinsic CurveCostProxy(X::ShimuraQuot, stage::MonStgElt) -> FldReElt
{A cheap upper-bound estimate of the cost of running `stage` on curve X, used only to
balance/order the parallel chunks.  Returns 0 for curves the stage skips.  The estimate
is the number of trace-formula terms a curve would need if it did NOT early-exit:
|W| * sum over good primes p<=Pbound, powers n=p^i, of sqrt(4*D*N*n).  It deliberately
over-estimates curves that early-exit (they then just finish fast), so the genuinely
heavy curves are always dispatched early.}
    R := RealField(6);
    if assigned X`IsSubhyp then return R!0; end if;
    g := X`g; DN := X`D * X`N; nW := #X`W;
    if g lt 3 then return R!0; end if;

    if stage in {"FilterByNonALInvolutions", "FilterByNonALInvolutionsStar"} then
        // cost ~ ModularSymbols(D*N); only curves with non-AL involutions do work, and only
        // those within the level cap (larger levels are skipped, so they cost ~nothing).
        if ((X`N mod 4 eq 0) or (Valuation(X`N, 3) eq 2)) and (DN le NonALModSymMaxLevel()) then
            return R!(DN*DN);
        end if;
        return R!0;
    end if;

    if stage eq "FilterByGeneralizedComplicatedFixedPoints" then
        // Only applicable curves (4|N or 9||N, within the level cap) do work; cost ~ the non-AL
        // fixed-point trace-formula sums over W.  Everything else the stage skips for free.
        if ((X`N mod 4 eq 0) or (Valuation(X`N, 3) eq 2)) and (DN le GeneralizedComplicatedMaxLevel()) then
            return R!(nW * g * Sqrt(R!DN));
        end if;
        return R!0;
    end if;

    if stage eq "FilterBySpecialFiber" then
        // D=1, D=6 and D=10 curves do work (a few primes, supersingular points, a Mobius check);
        // other discriminants are skipped for free.  Cost is roughly per-prime.  For D=6 the
        // supersingular points come from a degree ~p/24 hypergeometric polynomial whose roots
        // are found over F_{p^2}, so the dominant prime p|N enters linearly rather than via Sqrt.
        // For D=10 the cost is dominated by the Brandt module of discriminant 10p (~ class number,
        // linear in p) plus the Heun eigenvalue solve; weight it like D=6, slightly heavier.
        if X`D eq 1  then return R!(#PrimeDivisors(X`N) * Sqrt(R!X`N)); end if;
        if X`D eq 6  then return R!(#PrimeDivisors(X`N) * R!X`N); end if;
        if X`D eq 10 then return R!(#PrimeDivisors(X`N) * R!X`N * 2); end if;
        return R!0;
    end if;

    if stage in {"FilterByWeilPolynomial", "FilterByWeilPolynomialStar"} then
        // Weil cost is dominated by class-number streaming: each good prime p streams the
        // tables to depth 4*Qmax*p^g.  Bound p by the database, the affordability budget,
        // and the per-genus ceiling, then sum those depths over the good primes.
        ceil := AssociativeArray();
        ceil[3]:=53; ceil[4]:=53; ceil[5]:=37; ceil[6]:=29; ceil[7]:=23; ceil[8]:=17;
        Qmax := Max(X`W);
        b := WeilClassNumberPrimeBound(Qmax, g);
        bb := WeilBudgetPrimeBound(Qmax, g);
        if bb lt b then b := bb; end if;
        if IsDefined(ceil, g) and ceil[g] lt b then b := ceil[g]; end if;
        s := R!0;
        for p in PrimesUpTo(b) do
            if DN mod p eq 0 then continue; end if;
            s +:= R!(4*Qmax*p^g);
        end for;
        return s;
    end if;

    // trace-formula stages: choose prime bound and largest Hecke index per stage
    if stage in {"FilterByTrace", "FilterByTraceStar"} then
        Pbnd := 4*g^2;          // CheckHeckeTrace uses primes p <= 4g^2, n = p^v <= 4g^2
        Nmax := 4*g^2;
        imax := g;
    else
        return R!(nW * g * Sqrt(R!DN));    // generic light-stage fallback
    end if;

    s := R!0;
    for p in PrimesUpTo(Pbnd) do
        if DN mod p eq 0 then continue; end if;
        n := p; i := 1;
        while i le imax and (Nmax eq 0 or n le Nmax) do
            s +:= Sqrt(R!(4*DN*n));
            n *:= p; i +:= 1;
        end while;
    end for;
    return nW * s;
end intrinsic;
