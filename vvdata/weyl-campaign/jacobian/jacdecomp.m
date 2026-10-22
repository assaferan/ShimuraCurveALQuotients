// THE Q-ISOGENY DECOMPOSITION OF THE JACOBIAN OF EACH OPEN QUOTIENT, from the trace formula only.
//
//   magma -b LIST:=../nelson/nold_open61.txt jacdecomp.m < /dev/null
//
// WHY THE ZETA NUMERATOR SETTLES THIS.  If J ~ prod A_i over Q with A_i simple of dimension d_i,
// then P_p = prod P_{p,i} with deg P_{p,i} = 2 d_i.  A single P_{p,i} can factor further at a
// special prime (supersingular reduction, extra endomorphisms mod p), so the irreducible
// factorisation of P_p over Q is always a REFINEMENT of the true decomposition, never coarser.
// Hence the COARSEST pattern seen over several primes is the decomposition, and every prime that
// achieves it is a confirmation; the count of such primes is reported, because one prime alone
// could be an accident in the other direction -- a factor that happens not to split.
//
// This is the CHEAP method and it is EVIDENCE, not proof: a row whose coarsest pattern is realised
// at few primes deserves the exact computation in oldform.m, which decomposes over newform orbits
// and agrees with this on all 34 rows.
//
// For a simple factor of dimension 2 the endomorphism algebra is a real quadratic field, read off
// as the real quadratic subfield of Q[T]/(quartic).
//
// SELF-CONTAINED ON PURPOSE.  The zeta numerator is built here rather than by calling
// ZetaNumeratorViaTrace, which lives on the branch of PR #130 and is on neither main nor the
// campaign branch -- a script that cannot run on its own branch is useless.  The construction
// below is the same one: the point counts over F_p up to F_(p^g) give the power sums of the
// Frobenius eigenvalues, Newton's identities give the first g coefficients, and the functional
// equation gives the rest.  When #130 lands this can become a one-line call.

AttachSpec("ShimuraQuotients.spec");
SetColumns(0);

// P(T) for X over F_p.  The count over F_(p^(g+1)) is computed as well and compared with what P
// predicts: that comparison is the test, where checking the functional equation would prove
// nothing, since the upper coefficients are built from it and a wrong genus passes at every prime.
function ZetaNum(X, p, g)
    PT := PolynomialRing(Integers());
    if g eq 0 then return PT ! 1; end if;
    s := [ Rationals() | p^k + 1 - ComputePointsViaTrace(X, p, k) : k in [1..g] ];
    e := [ Rationals() | ];
    for k in [1..g] do
        acc := Rationals() ! 0;
        for j in [0..k-1] do
            ej := (j eq 0) select Rationals() ! 1 else e[j];
            acc +:= (-1)^(k-j-1) * ej * s[k-j];
        end for;
        Append(~e, acc/k);
    end for;
    cf := [ Rationals() | 1 ] cat [ (-1)^k * e[k] : k in [1..g] ];
    for k in [g+1..2*g] do Append(~cf, p^(k-g) * cf[2*g-k+1]); end for;
    error if not &and[ IsIntegral(c) : c in cf ],
        Sprintf("the counts at p = %o do not come from a curve of genus %o: zeta coefficients %o",
                p, g, cf);
    P := PT ! [ Integers() | c : c in cf ];
    ee := [ Rationals() | (-1)^k * Coefficient(P, k) : k in [1..2*g] ];
    ss := [ Rationals() | ];
    for k in [1..g+1] do
        acc := (-1)^(k-1) * k * ee[k];
        for i in [1..k-1] do acc +:= (-1)^(i-1) * ee[i] * ss[k-i]; end for;
        Append(~ss, acc);
    end for;
    error if ss[g+1] ne p^(g+1) + 1 - ComputePointsViaTrace(X, p, g+1),
        Sprintf("the numerator at p = %o built from counts up to F_(p^%o) does not predict the "
                cat "count over F_(p^%o)", p, g, g+1);
    return P;
end function;

lines := Split(Read(LIST), "\n");
printf "%-10o %-3o %-6o %-14o %-7o %o\n", "base", "g", "D*N", "decomposition", "primes",
       "fields of the dim-2 pieces";
for ln in lines do
    if #ln eq 0 or ln[1] eq "#" or Position(ln, "OPEN") eq 0 then continue; end if;
    fs := Split(ln, " ");
    base := fs[1];
    ps := Split(base, "_");
    D := StringToInteger(ps[1]); N := StringToInteger(ps[2]);
    i1 := Position(ln, "["); i2 := Position(ln, "]");
    ws := { StringToInteger(x) : x in Split(Substring(ln, i1+1, i2-i1-1), ",") };
    X := CreateShimuraQuot(D, N, ws);
    ok := true;
    try
        X`g := GenusShimuraCurveQuotient(D, N, ws);
    catch e
        ok := false;
    end try;
    if not ok then printf "%-10o  genus refused\n", base; continue; end if;
    X`CurveID := 0;
    pats := []; polys := [];
    for p in PrimesInInterval(3, 30) do
        if (D*N) mod p ne 0 then
            good := true; P := 0;
            try P := ZetaNum(X, p, X`g); catch e good := false; end try;
            if good then
                fa := Factorisation(P);
                dg := Sort(&cat[ [ Degree(t[1]) : i in [1..t[2]] ] : t in fa ]);
                Append(~pats, dg); Append(~polys, <p, fa, dg>);
            end if;
        end if;
    end for;
    if IsEmpty(pats) then printf "%-10o  no usable prime\n", base; continue; end if;
    best := pats[1];
    for q in pats do if #q lt #best then best := q; end if; end for;
    nhit := #[ q : q in pats | q eq best ];
    dims := [ d div 2 : d in best ];
    flds := {};
    for t in polys do
        if t[3] ne best then continue; end if;
        for u in t[2] do
            if Degree(u[1]) eq 4 then
                for sf in Subfields(NumberField(u[1]), 2) do
                    dd := Discriminant(MaximalOrder(sf[1]));
                    if dd gt 0 then Include(~flds, SquarefreeFactorization(dd)); end if;
                end for;
            end if;
        end for;
        break;
    end for;
    printf "%-10o %-3o %-6o %-14o %o/%-5o %o\n", base, X`g, D*N,
           Sprint(dims), nhit, #pats, IsEmpty(flds) select "-" else Sprint(flds);
end for;
quit;
