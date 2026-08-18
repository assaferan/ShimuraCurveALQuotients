// Why do the two 15_2 hauptmoduls need different m=0 multipliers (fs[-1] -> 4 validated,
// fs[-2] -> 2 required, code gives 4 for both)?
// Split data says fs[-1]'s 4 is ALL oo-block and fs[-2]'s 4 is ALL 0-block, so the 0-block
// assembly has never been validated by anything. This probe dumps the 0-block term by term:
// for every j with a nonzero cusp-0 coefficient, the bucket, its size, and each eta's contribution.
// Hypothesis under test: the bucket has size 2 (eta and -eta) and each gets full weight -> 2x.

AttachSpec("ShimuraQuotients.spec");
fh := Open("zeroblock_out.txt", "w");
emit := procedure(s) fprintf fh, "%o\n", s; Flush(fh); end procedure;

detail := procedure(foo, f0, Ldata, D, N, ~dummy)
    Q := Ldata`Q; disc_grp := Ldata`disc_grp; to_disc := Ldata`to_disc; denom := Ldata`denom;
    M := IsOdd(D*N) select 4*D*N else 2*D*N;
    Qint := ChangeRing(Q, Integers()); Qr := ChangeRing(Q, Rationals());
    negQ := -Qint;
    detprimes := Set(PrimeDivisors(Determinant(Qint)));
    Lfull := RSpaceWithBasis(IdentityMatrix(Integers(), 3));

    mod_M_to_vecs := AssociativeArray([0..M-1]);
    for j in [0..M-1] do mod_M_to_vecs[j] := []; end for;
    i0 := 0;
    for eta in disc_grp do
        if IsZero(eta) then i0 := eta; end if;
        v := ChangeRing(eta@@to_disc, Rationals());
        res := M*((v*Qr, v)/(2*denom^2));
        if not IsIntegral(res) then continue; end if;
        Append(~mod_M_to_vecs[Integers()!res mod M], eta);
    end for;

    contrib := function(eta, r, c)
        w_eta := ChangeRing(eta@@to_disc, Rationals())/denom;
        D0 := -(Numerator(r)*Denominator(r));
        K := QuadraticField(D0); dd := Discriminant(Integers(K));
        chi := KroneckerCharacter(dd);
        h := ClassNumber(K); wr := #TorsionSubgroup(UnitGroup(K));
        is_sq, cond_half := IsSquare(Rationals()!(r/AbsoluteValue(dd))); assert is_sq;
        Sc := Sort([p : p in detprimes join Set(PrimeDivisors(Numerator(r)))]);
        g := Rationals()!1;
        for p in Sc do
            g *:= LocalWhittakerAtOne(r, p, Vector(Rationals(), Eltseq(w_eta)), Lfull, negQ);
            if g eq 0 then break; end if;
        end for;
        if g eq 0 then return Rationals()!0; end if;
        en := &*[Rationals() | 1 - Evaluate(chi, p)/p : p in Sc];
        ed := &*[Rationals() | 1 - 1/(Rationals()!p)^2 : p in Sc];
        return c * cond_half * (Rationals()!h/wr) * (en/ed) * g;
    end function;

    pre := -96/(Rationals()!(D*N));

    emit(Sprintf("    M = %o, |disc_grp| = %o", M, #disc_grp));
    emit(Sprintf("    val(foo) = %o, val(f0) = %o", Valuation(foo), Valuation(f0)));
    // ---- oo-block, term by term
    emit("    -- oo-block (eta = 0) --");
    Too := Rationals()!0;
    for mm in [1..-Valuation(foo)] do
        c := Coefficient(foo, -mm);
        if c eq 0 then continue; end if;
        t := contrib(i0, Rationals()!mm, Rationals()!c);
        Too +:= t;
        emit(Sprintf("      m=%-4o c=%-8o term=%-10o  (x pre = %o)", mm, c, t, pre*t));
    end for;
    emit(Sprintf("    oo total = %o", pre*Too));

    // ---- 0-block, term by term, with bucket census
    emit("    -- 0-block --");
    T0 := Rationals()!0;
    for j in [1..-Valuation(f0)] do
        c := Coefficient(f0, -j);
        if c eq 0 then continue; end if;
        r := (Rationals()!j)/M;
        bucket := mod_M_to_vecs[j mod M];
        emit(Sprintf("      j=%-4o c=%-8o r=%-8o |bucket(%o)| = %o", j, c, r, j mod M, #bucket));
        sub := Rationals()!0;
        nz := 0;
        for eta in bucket do
            t := contrib(eta, r, Rationals()!c);
            sub +:= t;
            if t ne 0 then
                nz +:= 1;
                w_eta := ChangeRing(eta@@to_disc, Rationals())/denom;
                // is -eta also in the bucket, and is eta 2-torsion?
                two_tors := IsZero(eta + eta);
                emit(Sprintf("         eta = %-24o nonzero term = %-10o (x pre = %-8o) 2-torsion=%o",
                             Eltseq(w_eta), t, pre*t, two_tors));
            end if;
        end for;
        emit(Sprintf("        -> j=%o subtotal = %o  (nonzero etas: %o of %o)", j, pre*sub, nz, #bucket));
        T0 +:= sub;
    end for;
    emit(Sprintf("    0-block total = %o", pre*T0));
    emit(Sprintf("    GRAND TOTAL = %o", pre*(Too+T0)));
end procedure;

emit("START");
D := 15; N := 2;
Xstar := CreateShimuraQuot(D, N, Set(Divisors(D*N)));
Xstar`g := GenusShimuraCurveQuotient(D, N, Xstar`W); Xstar`CurveID := 0;
curves := GetQuotientsAndGenera([Xstar]);
_ := exists(star){c : c in curves | IsStarCurve(c)};
Ldata := ShimuraCurveLattice(D, N);
fs := BorcherdsForms(star, curves : Prec := 100);
emit(Sprintf("keys = %o", Sort([k : k in Keys(fs)])));
dummy := 0;
for k in [-1, -2] do
    emit(Sprintf("=========== form[%o] ===========", k));
    detail(qExpansionAtoo(fs[k], 1), qExpansionAt0(fs[k], 1), Ldata, D, N, ~dummy);
end for;
emit("DONE");
exit;
