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
        // cost ~ ModularSymbols(D*N); only curves with non-AL involutions do work
        if (X`N mod 4 eq 0) or (Valuation(X`N, 3) eq 2) then return R!(DN*DN); end if;
        return R!0;
    end if;

    // trace-formula stages: choose prime bound and largest Hecke index per stage
    if stage in {"FilterByTrace", "FilterByTraceStar"} then
        Pbnd := 4*g^2;          // CheckHeckeTrace uses primes p <= 4g^2, n = p^v <= 4g^2
        Nmax := 4*g^2;
        imax := g;
    elif stage eq "FilterByWeilPolynomial" then
        // primes up to the genus-scaled ceiling capped by the class-group-database bound
        ceil := AssociativeArray();
        ceil[3]:=53; ceil[4]:=53; ceil[5]:=37; ceil[6]:=29; ceil[7]:=23; ceil[8]:=17;
        b := WeilClassNumberPrimeBound(Max(X`W), g);
        if IsDefined(ceil, g) and ceil[g] lt b then b := ceil[g]; end if;
        Pbnd := b; Nmax := 0; imax := g;   // Nmax=0 -> no p^v<=Nmax cap, use i in [1..g]
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
