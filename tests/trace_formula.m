
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

procedure checkTracegDNew(g, D, N, k)
    assert GCD(D,N) eq 1;
    M := ModularSymbols(D*N,k,0);
    SDN := CuspidalSubspace(M);
    // We want the D-new subspace
    ps := PrimeDivisors(D);
    SDN_new := SDN;
    for p in ps do
        SDN_new := NewSubspace(SDN_new, p);
    end for;
    g_M := ActionOnModularSymbolsBasis(g,M);
    B := Matrix(Basis(VectorSpace(SDN_new)));
    trace := Trace(Solution(B, B*g_M));
    assert IsEven(Integers()!trace);
    from_formula := TraceFormulaGamma0gDNew(g,D,N,k);
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
    M2Q := MatrixAlgebra(Rationals(),2);
    M2Z := MatrixAlgebra(Integers(), 2);
    S2 := M2Q![2,1,0,2];
    return M2Z!(S2*W*S2^(-1));
end function;

function get_V3(N)
    alpha := Valuation(N,3);
    assert alpha eq 2;
    Q := 3^alpha;
    W := al_matrix(Q, N);
    M2Q := MatrixAlgebra(Rationals(),2);
    M2Z := MatrixAlgebra(Integers(), 2);
    S3 := M2Q![3,1,0,3];
    return M2Z!(S3*W*S3^(-1));
end function;

procedure testS2(bound, Dbound)
    Ns := [N : N in [1..bound] | N mod 8 eq 4];
    ks := [2,4,6];
    S2 := [2,1,0,2];
    for N in Ns do
        Ds := [D : D in [1..Dbound] | (MoebiusMu(D) eq 1) and (GCD(D, N) eq 1)];
        for k in ks do
            checkTraceg(S2, N, k);
            for D in Ds do
                checkTracegDNew(S2, D, N, k);
            end for;
        end for;
    end for;
end procedure;

procedure testV2(bound, Dbound)
    Ns := [N : N in [1..bound] | N mod 8 eq 0];
    ks := [2,4,6];
    for N in Ns do
        V2 := Eltseq(get_V2(N));
        Ds := [D : D in [1..Dbound] | (MoebiusMu(D) eq 1) and (GCD(D, N) eq 1)];
        for k in ks do
            checkTraceg(V2, N, k);
            for D in Ds do
                V2 := Eltseq(get_V2(D*N));
                checkTracegDNew(V2, D, N, k);
            end for;
        end for;
    end for;
end procedure;

procedure testV3(bound, Dbound)
    Ns := [N : N in [1..bound] | Valuation(N,3) eq 2];
    ks := [2, 4, 6];
    for N in Ns do
        V3 := Eltseq(get_V3(N));
        Ds := [D : D in [1..Dbound] | (MoebiusMu(D) eq 1) and (GCD(D, N) eq 1)];
        for k in ks do
            checkTraceg(V3, N, k);
            for D in Ds do
                V3 := Eltseq(get_V3(D*N));
                checkTracegDNew(V3, D, N, k);
            end for;
        end for;
    end for;
end procedure;