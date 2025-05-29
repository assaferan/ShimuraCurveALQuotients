
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

function TracegDNewModSym(g, D, N, k)
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
    return Integers()!(trace / 2);
end function;

procedure checkTracegDNew(g, g_subspaces, Q, D, N, k)
    trace := TracegDNewModSym(g, D, N, k);
    from_formula := TraceFormulaGamma0gDNew(g_subspaces,Q,D,N,k);
    // print "trace = ", trace, "from_formula = ", from_formula;
    assert trace eq from_formula;
end procedure;

procedure checkTraceVWDNew(p, Q, D, N, k)
    assert GCD(D,N) eq 1;
    M := ModularSymbols(D*N,k,0);
    SDN := CuspidalSubspace(M);
    // We want the D-new subspace
    ps := PrimeDivisors(D);
    SDN_new := SDN;
    for p in ps do
        SDN_new := NewSubspace(SDN_new, p);
    end for;
    if p eq 2 then
        Vp := get_V2(D*N);
    elif p eq 3 then
        Vp := get_V3(D*N);
    end if;
    W := al_matrix(Q, D*N);
    g := Eltseq(Vp*W);
    g_M := ActionOnModularSymbolsBasis(g,M);
    B := Matrix(Basis(VectorSpace(SDN_new)));
    trace := Trace(Solution(B, B*g_M));
    assert IsEven(Integers()!trace);
    from_formula := TraceFormulaGamma0VWDNew(p,Q,D,N,k);
    assert trace eq 2*from_formula;
end procedure;


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

function is_Hall_divisor(N, d)
    return Gcd(d, ExactQuotient(N, d)) eq 1;
end function;

function HallDivisors(N) 
    return [q : q in Divisors(N) | is_Hall_divisor(N, q)];
end function;

procedure testV2VW(bound, Dbound)
    Ns := [N : N in [1..bound] | N mod 8 eq 0];
    ks := [2,4,6];
    for N in Ns do
        Ds := [D : D in [1..Dbound] | (MoebiusMu(D) eq 1) and (GCD(D, N) eq 1)];
        for k in ks do
            for D in Ds do
                V2 := Eltseq(get_V2(D*N));
                hds := Set(HallDivisors(N));
                hds := {h : h in hds | GCD(h, 2) eq 1};
                Q := Random(hds);
                print "Q = ", Q;
                checkTraceVWDNew(2, Q, D, N, k);
            end for;
        end for;
    end for;
end procedure;

procedure testV3VW(bound, Dbound)
    Ns := [N : N in [1..bound] | Valuation(N,3) eq 2];
    ks := [2, 4, 6];
    for N in Ns do
        Ds := [D : D in [1..Dbound] | (MoebiusMu(D) eq 1) and (GCD(D, N) eq 1)];
        for k in ks do
            for D in Ds do
                hds := Set(HallDivisors(N));
                hds := {x : x in hds | x mod 3 eq 1};
                Q := Random(hds);
                print "Q = ", Q;
                checkTraceVWDNew(3, Q, D, N, k);
            end for;
        end for;
    end for;
end procedure;

// Testing Furumoto-Hasegawa theorem 4
procedure testFH4(curves : check_modular := CheckModularNonALInvolutionTrace)
    NWs := [
        < 63, {1,9}>,
        < 72, {1,9}>,
        <104, {1,104}>,
        <120, {1,5,24,120}>,
        <126, {1,9,7,63}>,
        <126, {1,9,14,126}>,
        <168, {1,24,56,21}>
    ];
    names := ["V3 W7", "V2 V3 W72", "V2 W13", "V2 W3", "V3 W1", "V3 W2", "V2 W3"];
    for idx in [1..#NWs] do
        NW := NWs[idx];
        assert exists(X){X : X in curves | X`D eq 1 and X`N eq NW[1] and X`W eq NW[2]};
        vname := names[idx];
        g, name, fix := check_modular(X);
        assert g eq 1;
        assert name eq vname; 
    end for;
    return;
end procedure;

function benchmark(curves, check_modular, ntimes)
    timings := [];
    Ns := Sort([N : N in {X`N : X in curves} | Valuation(N,2) ge 2 or Valuation(N,3) eq 2]);
    Ds := Sort([D : D in {X`D : X in curves}]);
    nNs := #Ns div 3;
    nDs := #Ds div 3;
    small_Ns := Ns[1..nNs];
    medium_Ns := Ns[nNs..2*nNs];
    big_Ns := Ns[2*nNs..#Ns];
    small_Ds := Ds[1..nDs];
    medium_Ds := Ds[nDs..2*nDs];
    big_Ds := Ds[2*nDs..#Ds];
    N_ranges := [small_Ns, medium_Ns, big_Ns];
    D_ranges := [small_Ds, medium_Ds, big_Ds];
    for D_range in D_ranges do
        time_by_D := [];
        for N_range in N_ranges do
            t := 0;
            DN_curves := [X : X in curves | X`N in N_range and X`D in D_range];
            if IsEmpty(DN_curves) then
                Append(~time_by_D, -1);
                continue;
            end if;
            for i in [1..ntimes] do
                X := Random(DN_curves);
                t0 := Cputime();
                b, name, fix := check_modular(X);
                t +:= Cputime() - t0;
            end for;
            Append(~time_by_D, t/ntimes);
        end for;
        Append(~timings, time_by_D);
    end for;
    return timings;
end function;