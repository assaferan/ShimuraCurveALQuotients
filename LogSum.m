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
            if (s`log_coeffs[p] eq 1) then coeff := ""; end if;
        end if;
        printf "%oLog%o", coeff, (p ne -1) select p else "oo";
    end for;
    return;
end intrinsic;

intrinsic RationalNumber(s::LogSm) -> FldRatElt
{.}
    require &and[IsIntegral(coeff) : coeff in s`log_coeffs] : "s does not represent a rational number!";
    if IsLogZero(s) then return 0; end if;
    if IsLogInfinity(s) then return Infinity(); end if;
    // Report WHICH prime has a runaway exponent. Magma's own failure here is
    // "Runtime error in '^': Argument 2 is too large", which names neither the prime nor the
    // exponent, and that is exactly how X_0^69(1) has been failing (it is the whole of that
    // base's triage record). The useful diagnostic is the offending (p, coeff) pair, not the
    // power that overflowed.
    //
    // ⇒ CAUSE FOUND 2026-09-14 (it is no longer "OPEN"; see HANDOFF.md, the night section, and run
    // with RUNAWAY=1 to reproduce the chain).  Three things have to line up:
    //   1. the Borcherds form has a huge principal part -- c(-m) reaches 19 digits at X_0^33(1)
    //      where X_0^21(1), which builds, tops out at 10 -- so every Schofer value carries a
    //      gigantic Log p component C;
    //   2. kappa_p(m) is IDENTICAL at every discriminant for exactly those m, so C is a pure
    //      COMMON factor that ought to cancel;
    //   3. ScaleForSchofer is NOT constant across the columns -- at d = -4 Ogg's condition halves
    //      W_size without n_d falling with it, giving scale -1/2 against -1/4 elsewhere -- so that
    //      column carries 2C where the rest carry C.
    // ReduceTable then subtracts the per-row MINIMUM, clearing every column but that one, which is
    // left holding the whole common factor.  C is an ARTIFACT: adding a kernel element changes the
    // form without changing its divisor, so C is not an invariant of the problem.  The fix to try
    // is to LLL-reduce the solution against the kernel to MINIMISE THE FORM'S coefficients -- note
    // the solution VECTOR is already tiny (maxsol 2 at 33_1); it is the echelon basis that is huge.
    //
    // ⚠ THE CAUSE IS NOT PRECISION, AND THIS COMMENT USED TO SAY IT WAS. "A runaway coefficient
    // means the Schofer sum diverged upstream" was an asserted cause, and it is REFUTED
    // (2026-09-13): re-running X_0^69(1) at Prec 300 instead of 100 reproduces the coefficient
    // 826241926712017437948244622352640031335552334419770916034895634120110682322 on Log23
    // BYTE-IDENTICALLY. A convergence or precision failure would move; an exact computation
    // returning a genuinely enormous integer does not. So this guard reports an OBSERVATION --
    // a coefficient past the threshold -- and the cause is open. Note the prime is RAMIFIED
    // (23 | 69), the same place condition 4's fractional exponents live.
    // ⇒ Do not spend a run raising Prec on this failure; that experiment is done.
    for p in Keys(s`log_coeffs) do
        error if AbsoluteValue(s`log_coeffs[p]) gt 10^5,
            Sprintf("RationalNumber: runaway log coefficient %o on Log%o (exceeds 10^5). "
                    * "Cause: huge principal part x a column whose ScaleForSchofer differs "
                    * "(see the note above; RUNAWAY=1 reproduces it). Full sum: %o",
                    s`log_coeffs[p], p, s);
    end for;
    ret := &*[Rationals() | p^(Integers()!s`log_coeffs[p]) : p in Keys(s`log_coeffs)];
    if (ret eq -1) then return Infinity(); end if;
    return ret;
end intrinsic;

intrinsic IsZero(s::LogSm) -> BoolElt
{.}
    reduce(s);
    return IsEmpty(Keys(s`log_coeffs));
end intrinsic;

intrinsic IsLogZero(s::LogSm) -> BoolElt
{.}
    return 0 in Keys(s`log_coeffs);
end intrinsic;

intrinsic IsLogInfinity(s::LogSm) -> BoolElt
{.}
    return -1 in Keys(s`log_coeffs);
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