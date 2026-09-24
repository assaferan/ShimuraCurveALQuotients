// DIRECT CHECK of the step my Schofer Thm 3.3 derivation turns on:
//   at m1 = 0 the conditions are  eta_+ + lambda_+ in L_+  AND  eta_- + lambda_- in L_-,
// and I claimed these force eta = 0 in L^v/L (because together they say eta + lambda lies in
// L_+ + L_- which sits inside L, while lambda is already in L).
// If that is right, only (eta, lambda) = (0, 0) survives.  If instead every ISOTROPIC eta keeps a
// lambda, the m = 0 term is a sum over the isotropic cosets -- which is what the oracle measures.
// Enumerate and count, using the pipeline's own Lplus/Lminus/L_quo construction.
AttachSpec("ShimuraQuotients.spec");

D := 15; N := 2;
if assigned DD then D := StringToInteger(DD); end if;
if assigned NN then N := StringToInteger(NN); end if;
Ld := ShimuraCurveLattice(D, N);
Q := ChangeRing(Ld`Q, Integers()); Qrat := ChangeRing(Q, Rationals());
dn := Ld`denom;

for d in [-7, -15, -40] do
    ok, lambda_v := true, 0;
    try
        lambda_v := ElementOfNorm(Ld`Q, -d, Ld`O, Ld`basis_L);
    catch e
        ok := false;
    end try;
    if not ok then printf "d = %o: no lambda\n", d; continue; end if;

    c_Lplus := Content(lambda_v);
    Lplus := RSpaceWithBasis(Matrix(lambda_v div c_Lplus));
    Lminus := Kernel(Transpose(Matrix(lambda_v*Q)));
    L := RSpaceWithBasis(IdentityMatrix(Integers(),3));
    L_quo, L_quo_map := L / (Lplus + Lminus);
    lambda_rat := ChangeRing(lambda_v, Rationals());
    nrm := (lambda_rat*Qrat, lambda_rat);

    // projections of a rational vector onto V_+ = Q.lambda and V_- = lambda^perp
    proj_plus := func<v | ((v*Qrat, lambda_rat)/nrm) * lambda_rat>;

    inLplus := func<v | IsCoercible(Lplus, [Rationals()!x : x in Eltseq(v)]) >;
    inLminus := func<v | IsCoercible(Lminus, [Rationals()!x : x in Eltseq(v)]) >;

    lams := [ l @@ L_quo_map : l in L_quo ];
    printf "\nd = %-5o |L/(L_+ + L_-)| = %-4o  |L^v/L| = %o\n", d, #L_quo, #Ld`disc_grp;
    cnt := AssociativeArray();
    total := 0; nonzero_eta := 0; iso_surviving := 0;
    for eta in Ld`disc_grp do
        v := ChangeRing(eta@@Ld`to_disc, Rationals())/dn;      // representative of eta in L^v
        nm := (v*Qrat, v)/2;
        c := 0;
        for lam in lams do
            w := v + ChangeRing(lam, Rationals());
            wp := proj_plus(w); wm := w - wp;
            if inLplus(wp) and inLminus(wm) then c +:= 1; end if;
        end for;
        if c gt 0 then
            total +:= c;
            if not IsZero(eta) then nonzero_eta +:= 1; end if;
            if IsIntegral(nm) then iso_surviving +:= 1; end if;
        end if;
        if not IsDefined(cnt, c) then cnt[c] := 0; end if;
        cnt[c] +:= 1;
    end for;
    printf "  surviving-lambda count per eta: %o\n",
           [<k, cnt[k]> : k in Sort([x : x in Keys(cnt)])];
    printf "  eta with at least one surviving lambda: total pairs %o, nonzero eta %o, isotropic %o\n",
           total, nonzero_eta, iso_surviving;
end for;
quit;
