// CONTROL for the incoherent hypothesis: how DISCRIMINATING is "|Diff_fin(r)| = 1"?
// If most r had |Diff| = 1 the rule would be near-vacuous, and its 8/9 score would mean nothing.
// The sharp question is whether it reproduces the LEVEL RULE: b(r) = 0 whenever N does not divide r.
AttachSpec("ShimuraQuotients.spec");

function is_square_Qp(a, p)
    v := Valuation(a, p);
    if IsOdd(v) then return false; end if;
    u := a / p^v; num := Numerator(u); den := Denominator(u);
    if p eq 2 then return ((num * InverseMod(den, 8)) mod 8) eq 1; end if;
    return IsSquare(GF(p) ! (num * InverseMod(den, p)));
end function;

function represents_Qp(diag3, a, p)
    diag := diag3 cat [-a];
    if not is_square_Qp(&*diag, p) then return true; end if;
    eps := &*[Integers() | HilbertSymbol(Rationals()!diag[i], Rationals()!diag[j], p)
              : i in [1..4], j in [1..4] | i lt j];
    return eps eq HilbertSymbol(Rationals()!(-1), Rationals()!(-1), p);
end function;

for b in [[15,2],[6,5],[10,3]] do
    D := b[1]; N := b[2];
    Ld := ShimuraCurveLattice(D, N);
    Dg := OrthogonalizeGram(ChangeRing(Ld`Q, Rationals()));
    diag3 := [Dg[i][i] : i in [1..3]];
    allp := Sort(Setseq({p : p in PrimeDivisors(2*3*5*7*11*13*17*19*23*Integers()!(&*diag3))}));
    printf "\n===== X0^%o(%o) =====\n", D, N;
    tab := AssociativeArray();
    for k in [0..4] do tab[k] := [0,0]; end for;   // [count with N|r, count with N not | r]
    ones_notdiv := [];
    for r in [1..60] do
        dset := [p : p in allp | not represents_Qp(diag3, Rationals()!r, p)];
        k := Minimum(#dset, 4);
        idx := (r mod N eq 0) select 1 else 2;
        tab[k][idx] +:= 1;
        if #dset eq 1 and (r mod N ne 0) then Append(~ones_notdiv, r); end if;
    end for;
    printf "  |Diff|   #(N | r)   #(N does not divide r)\n";
    for k in [0..4] do
        printf "  %-8o %-10o %o\n", k, tab[k][1], tab[k][2];
    end for;
    printf "  r with |Diff| = 1 but N does not divide r (these WOULD violate the level rule): %o\n",
           #ones_notdiv gt 12 select [#ones_notdiv] else ones_notdiv;
end for;
quit;
