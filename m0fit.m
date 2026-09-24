// Can the m = 0 multiplier be produced CHEAPLY, i.e. without the numeric coset sum?
//
// The oracle says  mult = (1/2) c_eta*(0)  at a nonzero isotropic eta*.  By Borcherds' obstruction,
// pairing F_f against a weight-3/2 holomorphic form G of the dual type gives
//     sum_{eta, m >= 0} c_eta(-m) b_eta(m) = 0,
// so taking G = the Eisenstein series whose constant term is e_eta* + e_-eta* (and using that all
// nonzero isotropic cosets carry the same constant term) gives
//     mult = (1/2) c_eta*(0) = -(1/4) sum_{eta, m>0} c_eta(-m) b^{eta*}_eta(m).
// m0_multiplier currently builds b from the Eisenstein series attached to the coset 0.  HYPOTHESIS:
// the local Whittaker should be evaluated at the coset SHIFTED by eta*.  Tested here as variants.
//
// The c_eta(-m) assignment itself (Lemma 24: the oo-block on eta = 0, the 0-block spread over the
// bucket) is NOT in question -- the numeric run verifies F_f's principal part to 1e-103.
AttachSpec("ShimuraQuotients.spec");

bases := [[15,2],[6,5],[10,3]];
// oracle multipliers, keys in the order Sort(Keys(BorcherdsForms)) = [-2,-1,9,10,11,12,13,14,15];
// entries flagged "false" were not independently measured by the sweep (they are predictions).
oracle := AssociativeArray();
oracle[[15,2]] := [<2,true>,<4,true>,<0,true>,<0,true>,<4,true>,<2,true>,<4,true>,<-2,true>,<2,true>];
oracle[[6,5]]  := [<0,false>,<0,false>,<6,true>,<3,true>,<12,false>,<9,false>,<3,true>,<0,true>,<3,true>];
oracle[[10,3]] := [<0,false>,<0,false>,<6,true>,<3,true>,<3,true>,<3,true>,<3,true>,<0,false>,<6,false>];

variants := ["zero", "minus", "plus"];

for b in bases do
    D := b[1]; N := b[2];
    Xstar := CreateShimuraQuot(D, N, Set(Divisors(D*N)));
    Xstar`g := GenusShimuraCurveQuotient(D, N, Xstar`W); Xstar`CurveID := 0;
    curves := GetQuotientsAndGenera([Xstar]);
    _ := exists(star){c : c in curves | IsStarCurve(c)};
    fs := BorcherdsForms(star, curves : Prec := 100);
    ks := Sort([k : k in Keys(fs)]);
    Ld := ShimuraCurveLattice(D, N);
    Q := Ld`Q; Qr := ChangeRing(Q, Rationals()); negQ := -ChangeRing(Q, Integers());
    dn := Ld`denom; M := IsOdd(D*N) select 4*D*N else 2*D*N;
    detprimes := Set(PrimeDivisors(Determinant(ChangeRing(Q, Integers()))));
    Lfull := RSpaceWithBasis(IdentityMatrix(Integers(), 3));

    // buckets, the trivial coset, and a nonzero isotropic coset eta*
    mod_M_to_vecs := AssociativeArray([0..M-1]);
    for j in [0..M-1] do mod_M_to_vecs[j] := []; end for;
    i0 := 0; isononzero := [];
    for eta in Ld`disc_grp do
        v := ChangeRing(eta@@Ld`to_disc, Rationals());
        nm := (v*Qr, v)/(2*dn^2);
        if IsZero(eta) then i0 := eta; end if;
        if IsIntegral(nm) and not IsZero(eta) then Append(~isononzero, eta); end if;
        if not IsIntegral(M*nm) then continue; end if;
        Append(~mod_M_to_vecs[Integers()!(M*nm) mod M], eta);
    end for;
    printf "\n======== X0^%o(%o)   #nonzero isotropic = %o ========\n", D, N, #isononzero;

    for vi->variant in variants do
      for star_idx in [1, 2] do
        if variant eq "zero" and star_idx eq 2 then continue; end if;
        wstar := ChangeRing(isononzero[star_idx]@@Ld`to_disc, Rationals())/dn;

        contrib := function(eta, r, c)
            w := ChangeRing(eta@@Ld`to_disc, Rationals())/dn;
            case variant:
                when "minus": w := w - wstar;
                when "plus":  w := w + wstar;
                else: ;
            end case;
            D0 := -(Numerator(r)*Denominator(r));
            K := QuadraticField(D0); dd := Discriminant(Integers(K));
            chi := KroneckerCharacter(dd);
            h := ClassNumber(K); wr := #TorsionSubgroup(UnitGroup(K));
            is_sq, cond_half := IsSquare(Rationals()!(r/AbsoluteValue(dd))); assert is_sq;
            Sc := Sort([p : p in detprimes join Set(PrimeDivisors(Numerator(r)))]);
            g := Rationals()!1;
            for p in Sc do
                g *:= LocalWhittakerAtOne(r, p, Vector(Rationals(), Eltseq(w)), Lfull, negQ);
            end for;
            en := &*[Rationals() | 1 - Evaluate(chi, p)/p : p in Sc];
            ed := &*[Rationals() | 1 - 1/(Rationals()!p)^2 : p in Sc];
            return c * cond_half * (Rationals()!h/wr) * (en/ed) * g;
        end function;

        Ts := [];
        for k in ks do
            foo := qExpansionAtoo(fs[k], 1); f0 := qExpansionAt0(fs[k], 1);
            T := Rationals()!0;
            for m in [1..-Valuation(foo)] do
                c := Coefficient(foo, -m);
                if c ne 0 then T +:= contrib(i0, Rationals()!m, Rationals()!c); end if;
            end for;
            for j in [1..-Valuation(f0)] do
                c := Coefficient(f0, -j);
                if c eq 0 then continue; end if;
                for eta in mod_M_to_vecs[j mod M] do
                    T +:= contrib(eta, (Rationals()!j)/M, Rationals()!c);
                end for;
            end for;
            Append(~Ts, T);
        end for;
        // does one constant per base send T to the oracle multiplier?
        ratios := [];
        ok := true;
        for i->k in ks do
            o := oracle[b][i][1];
            if Ts[i] eq 0 then
                if o ne 0 then ok := false; end if;
            else
                Append(~ratios, (Rationals()!o)/Ts[i]);
            end if;
        end for;
        consistent := ok and (#ratios eq 0 or #Set(ratios) eq 1);
        // also record the zero pattern, which is the sharpest structural signal
        zp := [Ts[i] eq 0 select 0 else 1 : i in [1..#ks]];
        zo := [oracle[b][i][1] eq 0 select 0 else 1 : i in [1..#ks]];
        printf "  %-6o star=%o  zeros %o vs oracle %o  %o  ratios=%o\n",
               variant, star_idx, zp, zo, (zp eq zo) select "ZEROS MATCH" else "zeros differ",
               Sprintf("%o distinct: %o", #Set(ratios),
                       #Set(ratios) le 4 select Setseq(Set(ratios)) else [Rationals()|]);
      end for;
    end for;
end for;
quit;
