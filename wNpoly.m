// W_N(s, m) for the N-SCALED binary lattice L_-, as a POLYNOMIAL in X = N^(-s).
// This is the raw material a derivation needs: order of vanishing at X = 1, the derivative there,
// and whether the polynomials fall into a recognisable closed form (the analogue for L_- of
// Kudla-Yang Prop 5.4 for the ternary Eichler lattice).
AttachSpec("ShimuraQuotients.spec");

for b in [<10,3,-8,[3,12,30]>, <15,2,-7,[2,10,30]>, <6,5,-19,[10,15,30]>] do
    D := b[1]; N := b[2]; d := b[3]; ms := b[4];
    Ld := ShimuraCurveLattice(D,N);
    Q := ChangeRing(Ld`Q, Integers()); Qrat := ChangeRing(Q, Rationals());
    lam := ElementOfNorm(Ld`Q, -d, Ld`O, Ld`basis_L);
    c_Lplus := Content(lam);
    Lminus := Kernel(Transpose(Matrix(lam*Q)));
    L := RSpaceWithBasis(IdentityMatrix(Integers(),3));
    Lplus := RSpaceWithBasis(Matrix(lam div c_Lplus));
    L_quo, L_quo_map := L / (Lplus + Lminus);
    lam_rat := ChangeRing(lam, Rationals());
    Bm := BasisMatrix(Lminus);
    Delta := Determinant(Bm*Q*Transpose(Bm));
    printf "\n===== X0^%o(%o) d=%o : disc(L_-) = %o, ord_%o = %o =====\n",
           D, N, d, Factorization(Delta), N, Valuation(Delta, N);
    seen := {};
    for m in ms do
        for mu_bar in L_quo do
            mu := mu_bar@@L_quo_map;
            c_mu := ((mu*Q, lam)/(lam*Q,lam));
            mu_minus := ChangeRing(mu, Rationals()) - c_mu*lam_rat;
            sqr_bd := m/(-d);
            for k in [Ceiling((-Sqrt(sqr_bd) - c_mu)*c_Lplus)..Floor((Sqrt(sqr_bd) - c_mu)*c_Lplus)] do
                x := (c_mu + k*c_Lplus^(-1))*lam_rat;
                mm := m - (x*Qrat,x)/2;
                if mm eq 0 then continue; end if;
                key := <mm, Eltseq(mu_minus)>;
                if key in seen then continue; end if;
                Include(~seen, key);
                W := LocalWhittakerPolynomial(mm, N, mu_minus, Lminus, Q);
                v1 := Evaluate(W, 1);
                dW := Derivative(W);
                ord := 0; WW := W;
                while Evaluate(WW, 1) eq 0 and ord lt 5 do WW := Derivative(WW); ord +:= 1; end while;
                printf "  m=%-4o mm=%-8o ord_%o(mm)=%-3o W_%o(X) = %-34o W(1)=%-8o W'(1)=%-8o ord_vanish=%o\n",
                       m, mm, N, Valuation(Numerator(mm),N) - Valuation(Denominator(mm),N),
                       N, W, v1, Evaluate(dW,1), ord;
            end for;
        end for;
    end for;
end for;
quit;
