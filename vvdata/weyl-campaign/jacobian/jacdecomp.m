// THE Q-ISOGENY DECOMPOSITION OF THE JACOBIAN OF EACH OPEN QUOTIENT, from the trace formula only.
//
//   magma -b LIST:=nold_open61.txt jacdecomp.m < /dev/null
//
// WHY THE ZETA NUMERATOR SETTLES THIS.  If J ~ prod A_i over Q with A_i simple of dimension d_i,
// then P_p = prod P_{p,i} with deg P_{p,i} = 2 d_i.  A single P_{p,i} can factor further at a
// special prime (supersingular reduction, extra endomorphisms mod p), so the irreducible
// factorisation of P_p over Q is always a REFINEMENT of the true decomposition, never coarser.
// Hence the COARSEST pattern seen over several primes is the decomposition, and every prime that
// achieves it is a confirmation; the count of such primes is reported, because one prime alone
// could be an accident in the other direction (a factor that happens not to split).
//
// For a simple factor of dimension 2 the endomorphism algebra is a real quadratic field, read off
// as the real quadratic subfield of Q[T]/(quartic) -- the same computation that gave the RM fields
// of the genus-two stars, which LMFDB confirmed 5/5.
AttachSpec("ShimuraQuotients.spec");
SetColumns(0);
lines := Split(Read(LIST), "\n");
printf "%-10o %-3o %-6o %-14o %-7o %o\n", "base", "g", "D*N", "decomposition", "primes", "fields of the dim-2 pieces";
for ln in lines do
    if #ln eq 0 or ln[1] eq "#" or Position(ln, "OPEN") eq 0 then continue; end if;
    fs := Split(ln, " ");
    base := fs[1];
    ps := Split(base, "_");
    D := StringToInteger(ps[1]); N := StringToInteger(ps[2]);
    // W is the bracketed list at the end of the row
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
        if (D*N) mod p eq 0 then continue; end if;
        good := true; P := 0;
        try P := ZetaNumeratorViaTrace(X, p); catch e good := false; end try;
        if not good then continue; end if;
        fa := Factorisation(P);
        dg := Sort(&cat[ [ Degree(t[1]) : i in [1..t[2]] ] : t in fa ]);
        Append(~pats, dg); Append(~polys, <p, fa, dg>);
    end for;
    if IsEmpty(pats) then printf "%-10o  no usable prime\n", base; continue; end if;
    // coarsest = fewest parts
    best := pats[1];
    for q in pats do if #q lt #best then best := q; end if; end for;
    nhit := #[ q : q in pats | q eq best ];
    dims := [ d div 2 : d in best ];
    // fields of the dimension-two pieces, taken at a prime realising the coarsest pattern
    flds := {};
    for t in polys do
        if t[3] ne best then continue; end if;
        for u in t[2] do
            if Degree(u[1]) ne 4 then continue; end if;
            for s in Subfields(NumberField(u[1]), 2) do
                dd := Discriminant(MaximalOrder(s[1]));
                if dd gt 0 then Include(~flds, SquarefreeFactorization(dd)); end if;
            end for;
        end for;
        break;
    end for;
    printf "%-10o %-3o %-6o %-14o %o/%-5o %o\n", base, X`g, D*N,
           Sprint(dims), nhit, #pats, IsEmpty(flds) select "-" else Sprint(flds);
end for;
quit;
