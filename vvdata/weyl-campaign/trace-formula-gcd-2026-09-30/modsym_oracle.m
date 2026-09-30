// Independent values for the TraceDNewALFixed pins, from modular symbols: the trace of T_n o w on
// the cuspidal subspace of level DN that is new at every p | D, averaged over w in W, with the
// Jacquet-Langlands sign (-1)^omega(gcd(w, D)).  Also genus(X_0(1848)/V2) via Tr(V2).
AttachSpec("ShimuraQuotients.spec");
import !"Geometry/ModSym/operators.m" : ActionOnModularSymbolsBasis;
function DNewTraceAvg(D, N, n, W)
    M := ModularSymbols(D*N, 2, 1);
    S := CuspidalSubspace(M);
    for p in PrimeDivisors(D) do S := NewSubspace(S, p); end for;
    if Dimension(S) eq 0 then return 0; end if;
    B := Matrix(Basis(VectorSpace(S)));
    T := HeckeOperator(M, n);
    tot := 0;
    for w in W do
        op := T * AtkinLehner(M, w);
        sgn := (-1)^#PrimeDivisors(GCD(w, D));
        tot +:= sgn * Trace(Solution(B, B*op));
    end for;
    return tot / #W;
end function;
W4 := {Integers() | 1, 8, 231, 1848};
W16 := {Integers() | 1, 6, 7, 10, 11, 15, 42, 66, 70, 77, 105, 110, 165, 462, 770, 1155};
Wa := {Integers() | 1, 5, 42, 210};
cases := [ <1, 1848, W4, [5, 25, 125, 7, 49]>, <6, 307, {Integers() | 1, 6}, [5, 25, 11]>,
           <10, 231, W16, [2, 4, 8, 13, 169]>, <14, 15, Wa, [11, 121, 1331]>,
           <1, 4522, {Integers() | 1, 4522}, [3, 9]>, <58, 7, {Integers() | 1, 14, 58, 203}, [3, 9, 27]> ];
for c in cases do
    D, N, W, ns := Explode(c);
    for n in ns do
        t0 := Cputime();
        v := DNewTraceAvg(D, N, n, W);
        f := TraceDNewALFixed(D, N, 2, n, W);
        printf "MODSYM <%o, %o, %o, W> = %o   formula %o   gcd(n,DN)=%o  %o   (%o s)\n", D, N, n, v, f, GCD(n, D*N),
            v eq f select "agree" else "DIFFER", RealField(3)!Cputime(t0);
    end for;
end for;
// genus of X_0(1848)/V2
M := ModularSymbols(1848, 2, 1); S := CuspidalSubspace(M);
g := Dimension(S);
V2 := get_V2(1848);
B := Matrix(Basis(VectorSpace(S)));
trV2 := Trace(Solution(B, B*ActionOnModularSymbolsBasis(Eltseq(V2), M)));
printf "GENUS X_0(1848): %o, Tr(V2) on S_2 = %o, genus of quotient = %o\n", g, trV2, (g + trV2)/2;
exit;
