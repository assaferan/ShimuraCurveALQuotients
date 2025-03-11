
import !"Geometry/ModSym/operators.m" : ActionOnModularSymbolsBasis;

procedure checkTraceg(g, N, k)
    M := ModularSymbols(N,k,0);
    C := CuspidalSubspace(M);
    g_M := ActionOnModularSymbolsBasis(g,M);
    B := Matrix(Basis(VectorSpace(C)));
    trace := Trace(Solution(B, B*g_M));
    assert IsEven(Integers()!trace);
    from_formula := TraceFormulaGamma0g(g,N,k);
    assert trace eq 2*from_formula;
end procedure;

function al_matrix(Q, N)
    d, x, y := XGCD(Q, N div Q);
    assert d eq 1;
    M2Z := MatrixAlgebra(Integers(), 2);
    return M2Z![Q*x, -y, N, Q];
end function;

function get_V2(N)
    alpha := Valuation(N,2);
    assert alpha ge 3;
    Q := 2^alpha;
    W := al_matrix(Q, N);
    M2Z := MatrixAlgebra(Integers(),2);
    S2 := M2Z![2,1,0,2];
    return S2*W*S2;
end function;

function get_V3(N)
    alpha := Valuation(N,3);
    assert alpha eq 2;
    Q := 3^alpha;
    W := al_matrix(Q, N);
    M2Z := MatrixAlgebra(Integers(),2);
    S3 := M2Z![3,1,0,3];
    return S3*W*S3^2;
end function;

procedure testS2(bound)
    Ns := [N : N in [1..bound] | N mod 8 eq 4];
    ks := [2,4,6];
    S2 := [2,1,0,2];
    for N in Ns do
        for k in ks do
            checkTraceg(S2, N, k);
        end for;
    end for;
end procedure;

procedure testV2(bound)
    Ns := [N : N in [1..bound] | N mod 8 eq 0];
    ks := [2,4,6];
    for N in Ns do
        V2 := Eltseq(get_V2(N));
        for k in ks do
            checkTraceg(V2, N, k);
        end for;
    end for;
end procedure;

procedure testV3(bound)
    Ns := [N : N in [1..bound] | Valuation(N,3) eq 2];
    ks := [2,4,6];
    for N in Ns do
        V3 := Eltseq(get_V3(N));
        for k in ks do
            checkTraceg(V3, N, k);
        end for;
    end for;
end procedure;