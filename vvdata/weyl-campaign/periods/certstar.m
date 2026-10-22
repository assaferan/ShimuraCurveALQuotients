// Certify a candidate star model against its newform: #C(F_p) = p + 1 - Tr(a_p).
// a_p is read from the modular symbols factor; taking it from the modular abelian variety instead
// is what made the in-line certification so costly that the level-2002 run was killed.
SetColumns(0);
PQ<x> := PolynomialRing(Rationals());
N := StringToInteger(LEVEL);
C := HyperellipticCurve(eval FPOL, eval HPOL);
printf "=== X_0^%o(1)*  candidate %o\n", N, C;
printf "genus %o ; discriminant supported on %o (D = %o)\n", Genus(C),
       PrimeDivisors(Integers() ! Discriminant(C)), PrimeDivisors(N);
MS := ModularSymbols(N, 2); NS := NewSubspace(CuspidalSubspace(MS));
B := 0; ps := PrimeDivisors(N);
for F in NewformDecomposition(NS) do
    if Dimension(F) ne 4 then continue; end if;
    if [ (AtkinLehnerOperator(F,q) eq -IdentityMatrix(Rationals(),4)) select -1 else 1 : q in ps ]
       eq [ -1 : q in ps ] then B := F; break; end if;
end for;
error if B cmpeq 0, "no star factor";
nm := 0; tot := 0; bad := [];
for p in PrimesInInterval(3, 60) do
    if N mod p eq 0 then continue; end if;
    ap := (Integers() ! Trace(HeckeOperator(B, p))) div 2;
    n := #Points(ChangeRing(C, GF(p)));
    tot +:= 1;
    if n eq p+1-ap then nm +:= 1; else Append(~bad, <p, n, p+1-ap>); end if;
end for;
printf "untwisted matches: %o of %o\n", nm, tot;
if nm eq tot then printf "CERTIFIED\n"; else printf "mismatches %o\n", bad; end if;
quit;
