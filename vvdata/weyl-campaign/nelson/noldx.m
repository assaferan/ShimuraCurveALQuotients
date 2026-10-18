// Extended count of the W-fixed differentials of X_0^D(N)/W, oldforms included.
//
// The D-new part of S_2(Gamma_0(DN)) is the direct sum, over M | N and over the newforms g of
// level exactly D*M, of the span of the pullbacks g | w_S, S a subset of the primes T dividing
// N/M.  On that span w_q for q | DM acts by the quaternionic eigenvalue eps_q(g) (a_g(q) at
// q | D, -a_g(q) at q | M) and w_p for p in T by the permutation S -> S xor {p}; so the
// eigenvectors are v_chi = sum_S chi(S) g|w_S for the characters chi of {+-1}^T, with w_m
// eigenvalue prod_{q | m, q | DM} eps_q(g) prod_{p | m, p in T} chi(p).  The W-fixed dimension
// is the number of (g, chi) with that product equal to 1 for every m in W, and must equal the
// genus of the quotient.
//
//   magma -b "BASES:=141_2:1,3,94,282;46_3:1,3" noldx.m < /dev/null
AttachSpec("/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/ShimuraQuotients.spec");
SetColumns(0);

// Per (D, N, W): newcnt, and a list of <M, label, chis> for the oldclasses that contribute.
function OldformCount(D, N, W)
    newcnt := 0; old := [* *];
    for M in Divisors(N) do
        T := PrimeDivisors(N div M);
        NFs := Newforms(CuspForms(Gamma0(D*M), 2));
        for orb in NFs do for f in orb do
            if Level(Parent(f)) ne D*M then continue; end if;
            qq := qExpansion(f, Max(PrimeDivisors(D*N)) + 2);
            chis := [];
            for bits in CartesianPower([1, -1], #T) do
                ok := true;
                for m in W do
                    if m eq 1 then continue; end if;
                    e := 1;
                    for p in PrimeDivisors(m) do
                        if D mod p eq 0 then e *:= Coefficient(qq, p);
                        elif M mod p eq 0 then e *:= -Coefficient(qq, p);
                        else e *:= bits[Index(T, p)]; end if;
                    end for;
                    if e ne 1 then ok := false; break; end if;
                end for;
                if ok then Append(~chis, [ bits[i] : i in [1..#T] ]); end if;
            end for;
            if M eq N then
                newcnt +:= #chis;
            elif #chis gt 0 then
                Append(~old, <M, [ Sprint(Coefficient(qq, n)) : n in [2..5] ], T, chis>);
            end if;
        end for; end for;
    end for;
    return newcnt, old;
end function;

for spec in Split(BASES, ";") do
    pp := Split(spec, ":");
    bn := Split(pp[1], "_");
    D := StringToInteger(bn[1]); N := StringToInteger(bn[2]);
    W := { StringToInteger(x) : x in Split(pp[2], ",") };
    g := GenusShimuraCurveQuotient(D, N, W);
    newcnt, old := OldformCount(D, N, W);
    oldcnt := &+([0] cat [ #o[4] : o in old ]);
    // does W pin the sign of any swap?  (a quotient where no m in W meets a prime of N/M for a
    // contributing oldclass tests nothing about the normalisation between the pullbacks)
    sens := exists{ o : o in old | exists{ m : m in W | exists{ p : p in o[3] | m mod p eq 0 } } };
    printf "%o_%o W=%o genus %o: new %o + old %o = %o  %o  swap-sign sensitive: %o\n",
        D, N, Sort(SetToSequence(W)), g, newcnt, oldcnt, newcnt + oldcnt,
        newcnt + oldcnt eq g select "MATCHES GENUS" else "*** MISMATCH ***", sens;
    for o in old do
        printf "      level %o (M=%o) a_2..a_5=%o  T=%o  chi=%o\n", D*o[1], o[1], o[2], o[3], o[4];
    end for;
    // genera of the quotients one step above: a 0 there makes the curve hyperelliptic for certain
    DN := D*N;
    hall := [ d : d in Divisors(DN) | d gt 1 and GCD(d, DN div d) eq 1 and d notin W ];
    above := [];
    for w in hall do
        W2 := W join { AtkinLehnerMul(w, m, DN) : m in W };
        Append(~above, <Sort(SetToSequence(W2)), GenusShimuraCurveQuotient(D, N, W2)>);
    end for;
    printf "      genera one step above: %o\n", [ a[2] : a in above ];
end for;
printf "NOLDX DONE\n";
quit;
