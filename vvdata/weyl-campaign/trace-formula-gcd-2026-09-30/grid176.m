// On the exact grid of tests/SfastUnitCondition.m blocks (2)/(3): how many full-space and newform
// traces does THIS tree get wrong vs modular symbols?  Plus TraceDNewALFixed(1,30,2,4,{1}).
AttachSpec("ShimuraQuotients.spec");
_ := ClassNumberLU(-4);
import "TraceFormula.m" : TraceFormulaGamma0HeckeAL, TraceFormulaGamma0HeckeALNew;
function ModSymTrace(N, k, n, Q : New := false)
    M := ModularSymbols(N, k, 1); C := CuspidalSubspace(M);
    if New then C := NewSubspace(C); end if;
    if Dimension(C) eq 0 then return 0; end if;
    B := Matrix(Basis(VectorSpace(C)));
    return Trace(Solution(B, B*(HeckeOperator(M, n) * AtkinLehner(M, Q))));
end function;
bf := 0; bn := 0; tot := 0;
for N in [12, 18, 20, 36, 45, 72] do
    for Q in [Q : Q in Divisors(N) | GCD(Q, N div Q) eq 1] do
        for n in [2, 3, 4, 5, 9] do
            if GCD(n, N) eq 1 then continue; end if;
            for k in [2, 4] do
                tot +:= 1;
                if TraceFormulaGamma0HeckeAL(N, k, n, Q) ne ModSymTrace(N, k, n, Q) then bf +:= 1; end if;
                if TraceFormulaGamma0HeckeALNew(N, k, n, Q) ne ModSymTrace(N, k, n, Q : New) then bn +:= 1; end if;
            end for;
        end for;
    end for;
end for;
printf "%o: on the test grid (%o tuples): full-space wrong %o, newform wrong %o; TraceDNewALFixed(1,30,2,4,{1}) = %o\n", tag, tot, bf, bn, TraceDNewALFixed(1, 30, 2, 4, {Integers() | 1});
exit;
