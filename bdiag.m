// Where does m0_multiplier's local assembly differ from the SOLVED Eisenstein coefficients?
// For each base and each index r in the level-rule support, print the exact b^{eta*}_0(r) recovered
// from the oracle, next to the code's per-unit weight
//     P(r, eta) = cond_half * (h/w) * (en/ed) * prod_{p in S_c} W_p(r, eta)
// at eta = 0 (what the code uses) and at a nonzero isotropic eta*, with the individual local factors
// broken out so it is visible WHICH prime is responsible.
AttachSpec("ShimuraQuotients.spec");

bases := [[15,2],[6,5],[10,3]];
btrue := AssociativeArray();
btrue[[15,2]] := [<2, 4>, <10, -4>, <30, 0>];
btrue[[6,5]]  := [<10, -6>, <15, -2>, <30, -6>];
btrue[[10,3]] := [<3, -2>, <12, -2>, <30, -6>];

for b in bases do
    D := b[1]; N := b[2];
    Ld := ShimuraCurveLattice(D, N);
    Q := Ld`Q; Qr := ChangeRing(Q, Rationals()); negQ := -ChangeRing(Q, Integers());
    dn := Ld`denom;
    detprimes := Set(PrimeDivisors(Determinant(ChangeRing(Q, Integers()))));
    Lfull := RSpaceWithBasis(IdentityMatrix(Integers(), 3));

    isononzero := [];
    for eta in Ld`disc_grp do
        v := ChangeRing(eta@@Ld`to_disc, Rationals());
        if IsIntegral((v*Qr, v)/(2*dn^2)) and not IsZero(eta) then Append(~isononzero, eta); end if;
    end for;
    zero := rep{eta : eta in Ld`disc_grp | IsZero(eta)};

    printf "\n============ X0^%o(%o)   D = %o, N = %o, det Q = %o ============\n",
           D, N, D, N, Determinant(ChangeRing(Q, Integers()));

    weight := function(r, eta)
        w := ChangeRing(eta@@Ld`to_disc, Rationals())/dn;
        D0 := -(Numerator(r)*Denominator(r));
        K := QuadraticField(D0); dd := Discriminant(Integers(K));
        chi := KroneckerCharacter(dd);
        h := ClassNumber(K); wr := #TorsionSubgroup(UnitGroup(K));
        is_sq, cond_half := IsSquare(Rationals()!(r/AbsoluteValue(dd))); assert is_sq;
        Sc := Sort([p : p in detprimes join Set(PrimeDivisors(Numerator(r)))]);
        locs := [];
        g := Rationals()!1;
        for p in Sc do
            lp := LocalWhittakerAtOne(r, p, Vector(Rationals(), Eltseq(w)), Lfull, negQ);
            Append(~locs, <p, lp>);
            g *:= lp;
        end for;
        en := &*[Rationals() | 1 - Evaluate(chi, p)/p : p in Sc];
        ed := &*[Rationals() | 1 - 1/(Rationals()!p)^2 : p in Sc];
        return cond_half * (Rationals()!h/wr) * (en/ed) * g, locs, dd, h, wr, cond_half;
    end function;

    for t in btrue[b] do
        r := Rationals()!t[1]; bt := t[2];
        P0, locs0, dd, h, wr, ch := weight(r, zero);
        Pst, locsst := weight(r, isononzero[1]);
        Pst2, _ := weight(r, isononzero[2]);
        printf "  r = %-4o b_true = %-5o | d = %-6o h/w = %-6o cond/2 = %-5o\n",
               t[1], bt, dd, Rationals()!h/wr, ch;
        printf "        P(r,0)    = %-12o locs = %o\n", P0, locs0;
        printf "        P(r,eta*) = %-12o locs = %o   [2nd eta*: %o]\n", Pst, locsst, Pst2;
        if P0 ne 0 then printf "        b_true/P(r,0)    = %o\n", bt/P0; end if;
        if Pst ne 0 then printf "        b_true/P(r,eta*) = %o\n", bt/Pst; end if;
    end for;
end for;
quit;
