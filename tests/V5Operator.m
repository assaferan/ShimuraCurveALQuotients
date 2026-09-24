// tests/V5Operator.m
//
// CI checks for the order-3 operator V_5 on X_0^D(N)/<w_25>, 25 || N (V5Operator.m).
//
// V_5 = S_5^{-1} T S_5 with T in Gamma_0(MD) of trace -1 mod 5 normalises Gamma_0(25MD)^+ but
// not Gamma_0(25MD), so it is computed on the modular symbols of Gamma_H(25MD),
// H = {d = +-1 mod 5}, and restricted to the Gamma_0, w_25 = +1 piece.  The checks below
// verify the structural claims that the (later) hyperellipticity test relies on:
//   [1] the Gamma_0 piece of S_2(Gamma_H)^{D-new} has dimension 2 g(X_0^D(N)),
//   [2] V_5^3 = 1, the w_25 = +1 piece is V_5-stable and the w_25 = -1 piece meets its image
//       trivially,
//   [3] on the w_25 = +1 piece, V_5 w_Q = w_Q V_5 when the 5-free part of Q is +-1 mod 5 and
//       V_5 w_Q = w_Q V_5^2 when it is +-2 mod 5 (these hold only modulo Gamma_0(N)^+, i.e.
//       only after restriction),
//   [4] tr V_5 = 0 on the complement of the w_25 = +1 piece,
//   [5] (dim + 2 tr V_5)/3 on the W-fixed piece is a non-negative even integer, so
//       g(X/<V_5>) is well defined,
//   [6] CanApplyV5 / V5OnQuotient agree with the genus formula for the quotient.

function hall_divisors(n)
    return [Q : Q in Divisors(n) | Q ne 1 and GCD(Q, n div Q) eq 1];
end function;

function five_free(Q)
    return Q div 5^Valuation(Q, 5);
end function;

function restrict_to(V, K)
    // matrix of the K-stable operator V (acting on row vectors) in a basis of K
    BK := BasisMatrix(K);
    return Matrix(Solution(BK, BK*V));
end function;

procedure testV5Matrix(D, M)
    DN := 25*M*D;
    V5 := get_V5(DN);
    assert Nrows(V5) eq 2 and BaseRing(V5) eq Integers();
    assert Determinant(V5) eq 25;
    // 5 V_5 = [5, i, -25MD, 5(1 - i M D)] with i M D = 3 mod 5
    i := V5[1,2];
    assert V5[1,1] eq 5 and V5[2,1] eq -DN and V5[2,2] eq 5*(1 - i*M*D);
    assert (i*M*D) mod 5 eq 3;
    // and T = S_5 V_5 S_5^{-1} is in Gamma_0(MD) with trace -1 mod 5
    M2Q := MatrixAlgebra(Rationals(), 2);
    S5 := M2Q![1, 0, 0, 1/5];
    T := S5 * (M2Q!V5/5) * S5^(-1);
    assert &and[IsIntegral(x) : x in Eltseq(T)] and Determinant(T) eq 1;
    assert Integers()!T[2,1] mod (M*D) eq 0 and Integers()!Trace(T) mod 5 eq 4;
end procedure;

procedure testV5Space(D, M)
    N := 25*M; DN := D*N;
    SH, B, K0 := GammaHSpaceV5(D, N);
    // [1] the Gamma_0 piece is the D-new cuspidal space of Gamma_0(DN)
    assert Dimension(K0) eq 2*GenusShimuraCurve(D, N);

    V := ActionOnSubspaceBasis(get_V5(DN), SH, B);
    assert V^3 eq 1;
    W25 := ActionOnSubspaceBasis(al_matrix(25, DN), SH, B);
    Kp := K0 meet Kernel(W25 - 1);
    Km := K0 meet Kernel(W25 + 1);
    assert Dimension(Kp) + Dimension(Km) eq Dimension(K0);
    // [2]
    assert Kp*V eq Kp;
    assert Dimension((Km*V) meet Km) eq 0;
    Vp := restrict_to(V, Kp);
    assert Vp^3 eq 1;
    // [4] trace zero on the complement
    assert Trace(V) eq Trace(Vp);
    // [5]
    fixdim := Integers()!(Dimension(Kp) + 2*Trace(Vp));
    assert fixdim ge 0 and fixdim mod 6 eq 0;
    // [3] relations with the Atkin-Lehner involutions on the w_25 = +1 piece
    for Q in hall_divisors(DN) do
        WQ := ActionOnSubspaceBasis(al_matrix(Q, DN), SH, B);
        assert Kp*WQ eq Kp;
        WQp := restrict_to(WQ, Kp);
        assert WQp^2 eq 1;
        if five_free(Q) mod 5 in {1, 4} then
            assert Vp*WQp eq WQp*Vp;
        else
            assert Vp*WQp eq WQp*Vp^2;
            // V_5 and w_Q then generate an S_3, so V_5 does not preserve the w_Q-eigenspaces
            // unless the S_3 acts through its abelianisation there.
            assert Vp*WQp ne WQp*Vp or Vp eq 1;
        end if;
    end for;
end procedure;

procedure testCanApplyV5()
    yes := [<1, 50, {1, 25}>, <1, 75, {1, 25}>, <1, 150, {1, 6, 25, 150}>, <6, 25, {1, 6, 25, 150}>,
            <1, 100, {1, 25}>, <1, 175, {1, 25}>];
    no := [<1, 50, {1, 50}>,            // w_25 not in W
           <1, 50, {1, 2, 25, 50}>,     // 2 = 2 mod 5
           <1, 150, {1, 2, 25, 50}>,
           <6, 25, {1, 2, 25, 50}>,
           <1, 150, {1, 3, 25, 75}>,    // 3 = -2 mod 5
           <1, 125, {1, 125}>,          // 25 does not exactly divide N
           <1, 30, {1, 5}>,             // 25 does not divide N
           <1, 50, {1}>];
    for t in yes do
        X := CreateShimuraQuot(t[1], t[2], t[3]);
        ok, _ := CanApplyV5(X);
        assert ok;
    end for;
    for t in no do
        X := CreateShimuraQuot(t[1], t[2], t[3]);
        ok, reason := CanApplyV5(X);
        assert not ok and Type(reason) eq MonStgElt and #reason gt 0;
    end for;
end procedure;

procedure testV5OnQuotient(D, N, W)
    X := CreateShimuraQuot(D, N, W);
    X`g := GenusShimuraCurveQuotient(D, N, W);
    V, K := V5OnQuotient(X);
    // [6] the W-fixed piece (with the D-sign convention) is 2 g(X) dimensional and V_5-stable
    assert Dimension(K) eq 2*X`g;
    assert Nrows(V) eq 2*X`g and V^3 eq 1;
    fixdim := Integers()!(2*X`g + 2*Trace(V));
    assert fixdim ge 0 and fixdim mod 6 eq 0;
end procedure;

procedure testV5OnQuotientRefuses(D, N, W)
    X := CreateShimuraQuot(D, N, W);
    X`g := GenusShimuraCurveQuotient(D, N, W);
    ok := false;
    try
        _ := V5OnQuotient(X);
    catch e
        ok := true;
    end try;
    assert ok;
end procedure;

// ── CI entry points (executed by run_tests.m) ───────────────────────────────
DMs := [<1, 2>, <1, 3>, <1, 4>, <1, 6>, <6, 1>, <1, 7>, <14, 1>, <1, 9>];
for DM in DMs do
    testV5Matrix(DM[1], DM[2]);
    testV5Space(DM[1], DM[2]);
end for;
testCanApplyV5();
for t in [<1, 50, {1, 25}>, <1, 75, {1, 25}>, <1, 100, {1, 25}>, <1, 150, {1, 25}>,
          <1, 150, {1, 6, 25, 150}>, <6, 25, {1, 25}>, <6, 25, {1, 6, 25, 150}>, <1, 175, {1, 25}>,
          <14, 25, {1, 14, 25, 350}>, <1, 225, {1, 9, 25, 225}>] do
    testV5OnQuotient(t[1], t[2], t[3]);
end for;
testV5OnQuotientRefuses(1, 50, {1, 50});
testV5OnQuotientRefuses(1, 150, {1, 2, 25, 50});
printf "V5Operator.m: all V_5 tests passed\n";
