// Step 1 of the beta != 0 implementation: enumerate the isotropic beta in A = L'/L,
// with their orders N_beta, and the data that controls how much of Schwagenscheidt Thm 1.4
// switches on at gamma = 0 (the block where the m=0 defect lives).
//
// At gamma = 0: (gamma,beta) = 0 so g = (N_beta, 0) = N_beta, hence N_g = N_beta, N_g' = 1,
//   eps_{beta,chi}(0) = G(chi)/N_beta   and   c^finite = prod_{p | N_beta} c^{finite,(p)}.
// Both are trivial exactly when N_beta = 1 (beta = 0) -- which is why every attempt so far
// collapsed to the same number. This probe reports the N_beta actually occurring, and how many
// Dirichlet characters mod N_beta of the right parity exist (and how many are primitive, since
// Thm 1.4 needs primitive chi and imprimitive ones need the Prop 1.3 oldform lifting).

AttachSpec("ShimuraQuotients.spec");
fh := Open("iso_out.txt", "w");
emit := procedure(s) fprintf fh, "%o\n", s; Flush(fh); end procedure;

emit("START");
for base in [<15,2>, <21,2>, <10,11>] do
    D, N := Explode(base);
    Ld := ShimuraCurveLattice(D, N);
    Q := Ld`Q; Qr := ChangeRing(Q, Rationals()); denom := Ld`denom;
    disc_grp := Ld`disc_grp; to_disc := Ld`to_disc;
    M := IsOdd(D*N) select 4*D*N else 2*D*N;
    emit(Sprintf("=== X0^%o(%o):  |A| = det Q = %o,  M = %o,  denom = %o ===",
                 D, N, Determinant(ChangeRing(Q,Integers())), M, denom));
    iso := [];
    nA := 0;
    for eta in disc_grp do
        nA +:= 1;
        v := ChangeRing(eta@@to_disc, Rationals());
        nm := (v*Qr, v)/(2*denom^2);
        if IsIntegral(nm) then
            w := v/denom;
            Nb := LCM([Denominator(x) : x in Eltseq(w)]);   // order of beta in L'/L
            Append(~iso, <eta, Nb, w>);
        end if;
    end for;
    emit(Sprintf("  |A| enumerated = %o,  #isotropic beta (Q(beta) in Z) = %o", nA, #iso));
    orders := {* t[2] : t in iso *};
    emit(Sprintf("  multiset of N_beta = %o", orders));
    for t in iso do
        eta, Nb, w := Explode(t);
        // characters mod N_beta: how many, how many primitive, and parity constraint
        nchi := EulerPhi(Nb);
        nprim := 0;
        if Nb eq 1 then
            nprim := 1;
        else
            G := DirichletGroup(Nb, CyclotomicField(Maximum(EulerPhi(Nb),1)));
            for c in Elements(G) do
                if Conductor(c) eq Nb then nprim +:= 1; end if;
            end for;
        end if;
        is2tors := IsZero(eta + eta);
        emit(Sprintf("    beta = %-22o  N_beta = %-4o  #chi mod N_beta = %-3o  #primitive = %-3o  2*beta=0? %o",
                     Eltseq(w), Nb, nchi, nprim, is2tors));
    end for;
end for;
emit("DONE");
