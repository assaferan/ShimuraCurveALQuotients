// Compute the m = 0 multiplier as Schofer's Theorem 3.3 actually gives it, and test its d-dependence.
//
// Theorem 3.3:  kappa_{gamma}(m) = sum_{mu_bar in L/(L_+ + L_-)} sum_{x in gamma_+ + mu_+ + L_+}
//                                     kappa^-_{gamma_- + mu_-}(m - Q(x)).
// At m = 0: Q(x) >= 0 forces Q(x) = 0, and since L_+ = L cap Q.lambda is RANK ONE POSITIVE DEFINITE
// (Q(x) = coeff^2 * (-d) with -d > 0), that means x = 0 exactly.  Combined with (11), i.e.
// kappa^-_mu(0) = 0 unless mu = 0, we get
//        kappa_eta(0) = kappa^-_0(0) * #{ mu_bar : 0 in eta_+ + mu_+ + L_+  and  eta_- + mu_- = 0 }
// so the multiplier of kappa^-_0(0) is   sum_eta c_eta(0) * (that count).
//
// The point of this probe: the splitting L_+ / L_- depends on lambda, hence on the CM DISCRIMINANT d.
// So this multiplier is d-DEPENDENT, whereas the shipped m0_multiplier is a single number per form.
//
// usage: magma -b DD:=<D> NN:=<N> m0sum.m

AttachSpec("ShimuraQuotients.spec");
SetVerbose("ShimuraQuotients", 0);

D := StringToInteger(DD); N := StringToInteger(NN);
Ld := ShimuraCurveLattice(D, N);
Q := Ld`Q; O := Ld`O; basis_L := Ld`basis_L;
disc_grp := Ld`disc_grp; to_disc := Ld`to_disc; denom := Ld`denom;
M := IsOdd(D*N) select 4*D*N else 2*D*N;
Qint := ChangeRing(Q, Integers()); Qrat := ChangeRing(Q, Rationals());

// buckets exactly as m0_multiplier / SchoferFormula0 form them
mod_M_to_vecs := AssociativeArray([0..M-1]);
for j in [0..M-1] do mod_M_to_vecs[j] := []; end for;
i0 := 0;
for eta in disc_grp do
    if IsZero(eta) then i0 := eta; end if;
    v := ChangeRing(eta@@to_disc, Rationals());
    res := M*((v*Qrat, v)/(2*denom^2));
    if not IsIntegral(res) then continue; end if;
    Append(~mod_M_to_vecs[Integers()!res mod M], eta);
end for;

// the count  #{ mu_bar : 0 in eta_+ + mu_+ + L_+ and eta_- + mu_- = 0 }  for one eta and one lambda
function eta_count(eta, lambda_v)
    c_Lplus := Content(lambda_v);
    Lplus := RSpaceWithBasis(Matrix(lambda_v div c_Lplus));
    Lminus := Kernel(Transpose(Matrix(lambda_v*Qint)));
    L := RSpaceWithBasis(IdentityMatrix(Integers(), 3));
    L_quo, L_quo_map := L / (Lplus + Lminus);
    lambda_rat := ChangeRing(lambda_v, Rationals());
    nrm := (lambda_rat*Qrat, lambda_rat);
    gamma := ChangeRing(eta@@to_disc, Rationals())/denom;
    c_gp := (gamma*Qrat, lambda_rat)/nrm;
    gamma_minus := gamma - c_gp*lambda_rat;
    cnt := 0;
    for mu_bar in L_quo do
        mu := ChangeRing(mu_bar@@L_quo_map, Rationals());
        c_mp := (mu*Qrat, lambda_rat)/nrm;
        mu_minus := mu - c_mp*lambda_rat;
        // x = 0 is attainable iff (c_gp + c_mp)*c_Lplus is an integer (then k realises it)
        if not IsIntegral((c_gp + c_mp)*c_Lplus) then continue; end if;
        if not IsZero(gamma_minus + mu_minus) then continue; end if;
        cnt +:= 1;
    end for;
    return cnt;
end function;

Xstar := CreateShimuraQuot(D, N, Set(Divisors(D*N)));
Xstar`g := GenusShimuraCurveQuotient(D, N, Xstar`W); Xstar`CurveID := 0;
curves := GetQuotientsAndGenera([Xstar]);
_ := exists(star){c : c in curves | IsStarCurve(c)};
fs := BorcherdsForms(star, curves : Prec := 100);
ks := Sort([k : k in Keys(fs)]);

// constant terms c_eta(0), in the pipeline's own convention
cconst := AssociativeArray();
for k in ks do
    foo := qExpansionAtoo(fs[k], 1); f0 := qExpansionAt0(fs[k], 1);
    ce := AssociativeArray();
    for eta in disc_grp do
        v := Rationals()!0;
        if eta eq i0 then v +:= Coefficient(foo, 0); end if;
        if eta in mod_M_to_vecs[0] then v +:= Coefficient(f0, 0); end if;
        ce[eta] := v;
    end for;
    cconst[k] := ce;
end for;

printf "[M0SUM] base D=%o N=%o  |disc_grp|=%o  M=%o\n", D, N, #disc_grp, M;
discs := [StringToInteger(x) : x in Split(DISCS, ",")];
for d in discs do
    lambda := ElementOfNorm(Q, -d, O, basis_L);
    cnts := AssociativeArray();
    for eta in disc_grp do cnts[eta] := eta_count(eta, lambda); end for;
    tot_supp := &+[Integers() | cnts[eta] : eta in disc_grp];
    for k in ks do
        mult := &+[Rationals() | cconst[k][eta] * cnts[eta] : eta in disc_grp];
        printf "[M0SUM] d=%-6o form=%-4o mult=%-10o (support=%o)\n", d, k, mult, tot_supp;
    end for;
end for;
exit;
