// Is the canonical representative CANONICAL as a FORM, or only at the cusps?
//
// Theorem 9.10 shows the CONSTANT TERMS of
//     E(D') := [ Theta(D',s) - Theta(D',1) ] / (mass_s - mass_1)
// are independent of D'.  That leaves open whether the FULL q-EXPANSIONS agree.
//   * if they do, W(m) = -a_E(m) is unambiguous -- a formula in representation
//     numbers of ternary quadratic forms, with no cusp-ideal caveat;
//   * if they differ, the difference is a CUSP form (equal constant terms), so
//     W(m) is defined only modulo the cusp ideal -- exactly the quotient the
//     paper describes -- and the functional is still well defined because cusp
//     forms pair to zero against principal parts (Bruinier-Funke).
// Either answer is worth having, and they are distinguishable in one run.
AttachSpec("ShimuraQuotients.spec");
DEPTH := 60; MMAX := DEPTH - 1;
R<q> := PowerSeriesRing(Rationals(), DEPTH);

grossTheta := function(Dp, Rl)   // mass-weighted AVERAGE theta, constant term 1
    Bq := QuaternionAlgebra(Dp); Oq := MaximalOrder(Bq);
    OR := Rl eq 1 select Oq else Order(Oq, Rl);
    CM := Matrix(Rationals(), 4, 4, &cat[ Eltseq(x) : x in Basis(OR) ]);
    LZ := Lattice(CM); Bv := [ Bq ! Eltseq(v) : v in Basis(LZ) ];
    den0 := Lcm([Denominator(Trace(x)) : x in Bv]);
    TrInt := Matrix(Integers(), #Bv, 1, [ Integers()!(den0*Trace(x)) : x in Bv ]);
    K := KernelMatrix(TrInt);
    S0 := [ &+[ K[i][j]*Bv[j] : j in [1..#Bv] ] : i in [1..Nrows(K)] ];
    GG := Matrix(Rationals(), 3,3, [ Trace(S0[i]*Conjugate(S0[j])) : j in [1..3], i in [1..3] ]);
    reps := Representatives(Genus(LatticeWithGram(GG)));
    num := [ Rationals()| 0 : m in [0..MMAX] ]; den := 0;
    for Lr in reps do
        wt := 1/#AutomorphismGroup(Lr);
        T := ThetaSeries(Lr, 2*MMAX+1);
        for m in [0..MMAX] do num[m+1] +:= wt*Coefficient(T, 2*m); end for;
        den +:= wt;
    end for;
    return [ x/den : x in num ], den;   // average theta, and the mass
end function;

mass := function(Dps, Rs)
    n := &*[ Integers() | p-1 : p in Dps ];
    m := Rationals()!n / (48*2^(#Dps-1));
    for x in Rs do m *:= Rationals()!(x+1)/2; end for;
    return m;
end function;

for D in [15, 39, 55] do
    pr := PrimeDivisors(D);
    printf "\n===== D = %o = %o, N = 1 =====\n", D, pr;
    E := [ ]; nm := [ ];
    for s in pr do
        Dps := [ p : p in pr | p ne s ]; Dp := &*Dps;
        m1 := mass(Dps, []); ms := mass(Dps, [s]);
        t1 := grossTheta(Dp, 1); ts := grossTheta(Dp, s);
        // [Theta_s - Theta_1]/(mass_s - mass_1) with Theta = mass * thetabar
        e := [ (ms*ts[i] - m1*t1[i])/(ms - m1) : i in [1..DEPTH] ];
        Append(~E, e); Append(~nm, Sprintf("(%o,1)+(%o,%o)", Dp, Dp, s));
        printf "  %-16o a(0..12) = %o\n", nm[#nm], [ e[i] : i in [1..13] ];
    end for;
    dif := [ E[1][i] - E[2][i] : i in [1..DEPTH] ];
    nz := [ i-1 : i in [1..DEPTH] | dif[i] ne 0 ];
    printf "  difference: %o\n", #nz eq 0 select "IDENTICALLY ZERO -- the representative is CANONICAL as a form"
        else Sprintf("nonzero at m = %o (first few coeffs %o) -- a CUSP form; W(m) lives in the quotient",
             [ nz[j] : j in [1..Minimum(8,#nz)] ], [ dif[nz[j]+1] : j in [1..Minimum(5,#nz)] ]);
end for;
printf "\nDONE\n";
quit;
