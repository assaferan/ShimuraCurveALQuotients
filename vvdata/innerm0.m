// Schofer: kappa_{eta,lambda}(m1) = sum_{x1 in eta_+ + lambda_+ + L_+} kappa_{eta_-+lambda_-}(m1 - Q(x1)),
// and kappa_{eta_-+lambda_-}(0) = k_0(0) * [eta_- + lambda_- in L_-].
// So whenever some x1 has Q(x1) = m EXACTLY, the m > 0 term picks up a genuine k_0(0) -- whose N-part
// is sum_{p | N/(N,d)} log p = log N at a FIRING disc.  That is an INNER m = 0, inside the x1-sum,
// not the outer m = 0 of sum_{m >= 0}.
// COUNT those (lambda, x1) for gamma = 0 (the Kappa0 case) and compare with the target A_m.
AttachSpec("ShimuraQuotients.spec");
D := 15; N := 2; d := -7;
if assigned DD then D := StringToInteger(DD); end if;
if assigned NN then N := StringToInteger(NN); end if;
if assigned DISC then d := StringToInteger(DISC); end if;
targets := AssociativeArray();
targets[[15,2]] := [<2,-1>,<10,1>,<30,0>,<1,0>,<3,0>,<5,0>,<15,0>];
targets[[6,5]]  := [<10,3/2>,<15,1/2>,<30,3/2>,<1,0>,<2,0>,<3,0>,<6,0>];
targets[[10,3]] := [<3,1/2>,<12,1/2>,<30,3/2>,<1,0>,<2,0>,<5,0>,<10,0>];

Ld := ShimuraCurveLattice(D,N);
Q := ChangeRing(Ld`Q, Integers()); Qrat := ChangeRing(Q, Rationals());
lambda_v := ElementOfNorm(Ld`Q, -d, Ld`O, Ld`basis_L);
c_Lplus := Content(lambda_v);
Lplus := RSpaceWithBasis(Matrix(lambda_v div c_Lplus));
Lminus := Kernel(Transpose(Matrix(lambda_v*Q)));
L := RSpaceWithBasis(IdentityMatrix(Integers(),3));
L_quo, L_quo_map := L / (Lplus + Lminus);
lambda_rat := ChangeRing(lambda_v, Rationals());
printf "X0^%o(%o) at d = %o : |L/(L_+ + L_-)| = %o, N-part of k_0(0) = %o * log %o\n",
       D, N, d, #L_quo,
       (FundamentalDiscriminant(d) mod N eq 0) select 0 else 1, N;
printf "  %-5o %-10o %-10o %o\n", "m", "count", "count/2", "target A_m";
for t in targets[[D,N]] do
    m := t[1]; want := t[2];
    cnt := 0;
    for mu_bar in L_quo do
        mu := mu_bar@@L_quo_map;
        c_mu_plus := ((mu*Q, lambda_v)/(lambda_v*Q,lambda_v));
        mu_plus := c_mu_plus*lambda_rat;
        mu_minus := ChangeRing(mu, Rationals()) - mu_plus;
        // kappa^-(0) is nonzero only when mu_minus lies in L_minus
        if not IsCoercible(Lminus, [Rationals()!x : x in Eltseq(mu_minus)]) then continue; end if;
        sqr_bd := m/(-d);
        lb := Ceiling((-Sqrt(sqr_bd) - c_mu_plus)*c_Lplus);
        ub := Floor((Sqrt(sqr_bd) - c_mu_plus)*c_Lplus);
        for k in [lb..ub] do
            x := (c_mu_plus + k*c_Lplus^(-1))*lambda_rat;
            if (m - (x*Qrat,x)/2) eq 0 then cnt +:= 1; end if;
        end for;
    end for;
    printf "  %-5o %-10o %-10o %o\n", m, cnt, cnt/2, want;
end for;
quit;
