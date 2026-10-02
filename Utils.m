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
{A cheap estimate of the cost of running `stage` on curve X, used only to order the parallel
chunks heavy-first.  Returns 0 for curves the stage skips.  For the Weil-polynomial stages it
is the sum of p^g over the good primes the stage will use; for the trace stages it is the
number of trace-formula terms a curve would need if it did not early-exit, |W| times the sum of
sqrt(4*D*N*n) over the primes and powers n = p^i; the other stages have their own rough
measures below.  It deliberately over-estimates curves that early-exit (they then just finish
fast), so the genuinely heavy curves are always dispatched early.}
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

    if stage eq "FilterByAutomorphismGroup" then
        // Cheap coset enumeration; the occasional expensive part is TraceDNewQuotient for a
        // non-AL central involution, which only arises with S2/V2/V3 (4 | N or 9 || N).
        if (X`N mod 4 eq 0) or (Valuation(X`N, 3) eq 2) then return R!(g * Sqrt(R!DN)); end if;
        return R!1;
    end if;

    if stage in {"FilterByTwistedTrace", "FilterByTwistedWeilPolynomial",
                 "FilterByTwistedTraceStar", "FilterByTwistedWeilPolynomialStar"} then
        // Modular symbols of level D*N plus Hecke operators: grows roughly like DN^2.  The worker
        // takes the maximum over a level (the modular symbols are shared by its curves).
        if stage in {"FilterByTwistedWeilPolynomial", "FilterByTwistedWeilPolynomialStar"} and
           not ((g eq 3) or (g eq 4 and exists{p : p in [2,3,5] | DN mod p ne 0}) or (g in {5,6} and IsOdd(DN))) then
            return R!0;
        end if;
        // no admissible involution h != 1 (always so for a star curve without V2/V3): skipped
        has_ops := (#{Q : Q in Divisors(DN) | GCD(Q, DN div Q) eq 1} gt nW)
                   or ((X`N mod 4 eq 0) and &and[IsOdd(w) : w in X`W])
                   or (X`N mod 8 eq 0) or ((Valuation(X`N, 3) eq 2) and (9 in X`W));
        if not has_ops then return R!0; end if;
        return R!(DN) * R!(DN);
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
        if X`D eq 22 then return R!(#PrimeDivisors(X`N) * R!X`N * 2); end if;   // Brandt(22p) + CM lift
        return R!0;
    end if;

    if stage in {"FilterByWeilPolynomial", "FilterByWeilPolynomialStar"} then
        // The stage asks for #X(F_{p^v}), v <= g, at every good prime p up to its bound (the
        // class-number database, the disc budget, and the per-genus ceiling).  The dominant call
        // is the trace of T_n at n = p^g, which Eichler-Selberg evaluates as a sum over the
        // elements w of W of a sum over t with t^2 < 4n/Q_w, i.e. about 2 sqrt(n/Q_w) terms, so
        // the term count per prime grows like p^(g/2) and the number of inner sums like #W.
        //
        // Measured on lovelace with the class-number tables, over 36 curves spanning genus 3 to
        // 7, #W from 1 to 64 and level from 30 to 30030, on main at 95b19e6 (so with #56, #57
        // and #58): vvdata/weyl-campaign/weil-retime-2026-10-02/ on the m0-theta-campaign
        // branch.  Of the 271 curve pairs whose times differ by more than a factor 10, this
        // orders 246 correctly against 198 for the sum of p^g that the function used before; of
        // the 520 pairs differing by more than a factor 2, 440 against 343.  So both the #W and
        // the exponent g/2 in place of g are improvements, and that is as much as the data
        // supports: among #W * sum p^(a g) for a between 0.5 and 0.8, and #W * max p^(a g), the
        // counts differ by less than the noise of 36 curves, so the term count is kept because
        // it is the principled one rather than the best-scoring one.
        //
        // ⚠ A residual spread of 237x remains (time over estimate), so this orders chunks and is
        // never a cost.  It still mis-orders real pairs: X_0^6(97)/W2 (g = 6, 74 min) sits below
        // X_0^210(73)/W32 (g = 4, 86 s).  The two slowest curves of the 36 are X_0(595)/W8
        // (g = 5, 101 min) and that X_0^6(97)/W2, both at level about 590 with a SMALL W, so the
        // heavy shape is high genus at moderate level and a large W is not what makes a curve
        // slow -- #W earns its place in the ordering, not in the extremes.  The level is weakly
        // positive (fitting #W^a (sum p^(g/2))^b (D*N)^c gives c = +0.15), not absent.
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
            s +:= Sqrt(R!(p^g));
        end for;
        return R!nW * s;
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
