// Validation of the generalized complicated-fixed-point filter.
//   Part A: NumFixedPointsNonALOnX (trace formula) vs an independent modular-symbols
//           computation of the same fixed-point count (D=1, where the D-new subspace is
//           the full cuspidal space and the Shimura sign twist is trivial).
//   Part B: run the criterion on real curves; it must NEVER claim "non-hyperelliptic" for
//           a tuple already proven HYPERELLIPTIC by another stage (no contradictions).
SetQuitOnError(true);
AttachSpec("ShimuraQuotients.spec");
Attach("GeneralizedComplicatedFixedPoints.m");
SetVerbose("ShimuraQuotients", 0);
import !"Geometry/ModSym/operators.m" : ActionOnModularSymbolsBasis;

// Trace of the matrix g on the D-new holomorphic differentials, via modular symbols.
// (The 2g-dimensional modular-symbols trace is 2t; this returns t.)
function ModSymTraceOnForms(g, D, N)
    M := ModularSymbols(D*N, 2, 0);
    SDN := CuspidalSubspace(M);
    SDN_new := SDN;
    for p in PrimeDivisors(D) do SDN_new := NewSubspace(SDN_new, p); end for;
    g_M := ActionOnModularSymbolsBasis(g, M);
    B := Matrix(Basis(VectorSpace(SDN_new)));
    trace := Trace(Solution(B, B*g_M));
    assert IsEven(Integers()!trace);
    return Integers()!(trace/2);
end function;

// #Fix of an involution = 2 - 2*trace(on differentials)   (Lefschetz / Riemann-Hurwitz).
function ModSymNuOnX(V, Q, D, N)
    t := ModSymTraceOnForms(Eltseq(V * al_matrix(Q, D*N)), D, N);
    return 2 - 2*t;
end function;

print "=== PART A: helper vs modular-symbols cross-check (D=1) ===";
S2 := Matrix(Integers(),2,2,[2,1,0,2]);
nA := 0;
for tup in [ <36,"S2">, <36,"V3">, <72,"S2">, <72,"V2">, <72,"V3">,
             <32,"V2">, <32,"S2">, <144,"V2">, <144,"V3">, <100,"S2"> ] do
    N := tup[1]; vname := tup[2]; DN := N;
    // availability guards: S2 needs 4|N, V2 needs 8|N, V3 needs 9||N
    if (vname eq "S2") and (N mod 4 ne 0) then continue; end if;
    if (vname eq "V2") and (N mod 8 ne 0) then continue; end if;
    if (vname eq "V3") and (Valuation(N,3) ne 2) then continue; end if;
    if   vname eq "S2" then V := S2;
    elif vname eq "V2" then V := get_V2(DN);
    else V := get_V3(DN); end if;
    for Q in [q : q in Divisors(DN) | GCD(q, DN div q) eq 1] do
        if (vname eq "S2") and IsEven(Q) then continue; end if;
        if (vname eq "V3") and
           exists(p){p : p in PrimeDivisors(Q) | (p^Valuation(Q,p) mod 3) eq 2} then continue; end if;
        h := NumFixedPointsNonALOnX(V, vname, Q, 1, N);
        m := ModSymNuOnX(V, Q, 1, N);
        if h ne m then
            printf "  *** MISMATCH N=%o %o W%o : helper=%o modsym=%o\n", N, vname, Q, h, m;
        end if;
        assert h eq m;
        nA +:= 1;
    end for;
end for;
printf "  %o cross-checks, ALL PASSED.\n", nA;

print "";
print "=== PART B: no-contradictions scan on real curves ===";
curves := eval Read("data/par/curves_after_UpdateCurves7.dat");
printf "loaded %o curves\n", #curves;

known := AssociativeArray(); rep := AssociativeArray(); keys := [];
for X in curves do
    key := <X`D, X`N, X`W>;
    if not IsDefined(rep, key) then rep[key] := X; known[key] := "open"; Append(~keys, key); end if;
    if assigned X`IsHyp and X`IsHyp then known[key] := "hyp"; end if;
    if (known[key] ne "hyp") and assigned X`IsSubhyp and (not X`IsSubhyp) then known[key] := "nonhyp"; end if;
end for;
cap := GeneralizedComplicatedMaxLevel();
applicable := [k : k in keys |
    ((k[2] mod 4 eq 0) or (Valuation(k[2],3) eq 2)) and (rep[k]`g ge 3) and (k[1]*k[2] le cap)];
printf "unique tuples: %o ; applicable (4|N or 9||N, g>=3, DN<=%o): %o\n", #keys, cap, #applicable;

flagged := 0; contradictions := 0; consistent := 0; newp := 0; examples := [];
for idx->k in applicable do
    ok, witness := CheckGeneralizedComplicatedFixedPoints(rep[k]);
    if ok then
        flagged +:= 1;
        if known[k] eq "hyp" then
            contradictions +:= 1;
            printf "  *** CONTRADICTION: D=%o N=%o g=%o known HYP but flagged: %o\n",
                k[1], k[2], rep[k]`g, witness;
        elif known[k] eq "nonhyp" then consistent +:= 1;
        else
            newp +:= 1;
            if #examples lt 15 then Append(~examples, Sprintf("D=%o N=%o g=%o :: %o", k[1], k[2], rep[k]`g, witness)); end if;
        end if;
    end if;
    if (idx mod 500 eq 0) then printf "  ...%o/%o\n", idx, #applicable; end if;
end for;

print "---------------- RESULT ----------------";
printf "flagged non-hyperelliptic    : %o\n", flagged;
printf "  CONTRADICTIONS (must be 0) : %o\n", contradictions;
printf "  consistent w/ known nonhyp : %o\n", consistent;
printf "  NEW prunings (were open)   : %o\n", newp;
print "---- example NEW prunings ----";
for e in examples do print "  ", e; end for;
quit;
