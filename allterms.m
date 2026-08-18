// Full per-term dump of the m=0 multiplier across ALL BorcherdsForms of a base, with per-prime
// local factors, so the 4-vs-2 discrepancy can be fit against every available ground truth:
//   15_2 fs[-1] must total 4 (Table-45 validated), fs[-2] must total 2 (hauptmodul relation),
//   10_11 every form must total 0 (main passes those tests with no m=0 term).
// Emits one line per (block, m/j, eta) contribution with its ingredients.

AttachSpec("ShimuraQuotients.spec");
outname := "allterms_out.txt";
if assigned OUTNAME then outname := OUTNAME; end if;
fh := Open(outname, "w");
emit := procedure(s) fprintf fh, "%o\n", s; Flush(fh); end procedure;

dump_form := procedure(foo, f0, Ldata, D, N)
    Q := Ldata`Q; disc_grp := Ldata`disc_grp; to_disc := Ldata`to_disc; denom := Ldata`denom;
    M := IsOdd(D*N) select 4*D*N else 2*D*N;
    Qint := ChangeRing(Q, Integers()); Qr := ChangeRing(Q, Rationals());
    negQ := -Qint;
    dQ := Determinant(Qint);
    detprimes := Set(PrimeDivisors(dQ));
    Lfull := RSpaceWithBasis(IdentityMatrix(Integers(), 3));
    pre := -96/(Rationals()!(D*N));

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

    // returns term and a descriptive string of its ingredients
    contrib := function(eta, r, c)
        w_eta := ChangeRing(eta@@to_disc, Rationals())/denom;
        D0 := -(Numerator(r)*Denominator(r));
        K := QuadraticField(D0); dd := Discriminant(Integers(K));
        chi := KroneckerCharacter(dd);
        h := ClassNumber(K); wr := #TorsionSubgroup(UnitGroup(K));
        is_sq, cond_half := IsSquare(Rationals()!(r/AbsoluteValue(dd))); assert is_sq;
        Sc := Sort([p : p in detprimes join Set(PrimeDivisors(Numerator(r)))]);
        g := Rationals()!1;
        locs := [];
        for p in Sc do
            lp := LocalWhittakerAtOne(r, p, Vector(Rationals(), Eltseq(w_eta)), Lfull, negQ);
            Append(~locs, <p, lp, Valuation(Numerator(r), p) - Valuation(Denominator(r), p)>);
            g *:= lp;
        end for;
        en := &*[Rationals() | 1 - Evaluate(chi, p)/p : p in Sc];
        ed := &*[Rationals() | 1 - 1/(Rationals()!p)^2 : p in Sc];
        t := c * cond_half * (Rationals()!h/wr) * (en/ed) * g;
        desc := Sprintf("d=%o h=%o w=%o cond/2=%o en/ed=%o g=%o locs(p,W,v_p(r))=%o",
                        dd, h, wr, cond_half, en/ed, g, locs);
        return t, desc;
    end function;

    Too := Rationals()!0; T0 := Rationals()!0;
    for mm in [1..-Valuation(foo)] do
        c := Coefficient(foo, -mm);
        if c eq 0 then continue; end if;
        t, desc := contrib(i0, Rationals()!mm, Rationals()!c);
        Too +:= t;
        emit(Sprintf("    oo m=%-5o c=%-6o contrib=%-8o | %o", mm, c, pre*t, desc));
    end for;
    for j in [1..-Valuation(f0)] do
        c := Coefficient(f0, -j);
        if c eq 0 then continue; end if;
        r := (Rationals()!j)/M;
        bucket := mod_M_to_vecs[j mod M];
        nzs := [];
        sub := Rationals()!0;
        for eta in bucket do
            t, desc := contrib(eta, r, Rationals()!c);
            sub +:= t;
            if t ne 0 then
                w_eta := ChangeRing(eta@@to_disc, Rationals())/denom;
                Append(~nzs, <Eltseq(w_eta), pre*t, IsZero(eta+eta), desc>);
            end if;
        end for;
        T0 +:= sub;
        emit(Sprintf("    0  j=%-5o c=%-6o r=%-7o |bkt|=%-4o subtotal=%-8o", j, c, r, #bucket, pre*sub));
        for z in nzs do
            emit(Sprintf("         eta=%-22o term=%-8o 2tors=%-6o | %o", z[1], z[2], z[3], z[4]));
        end for;
    end for;
    emit(Sprintf("    ==> oo = %o , zero = %o , TOTAL = %o", pre*Too, pre*T0, pre*(Too+T0)));
end procedure;

emit("START");
D := StringToInteger(DD); N := StringToInteger(NN);
emit(Sprintf("=========== X0^%o(%o) ===========", D, N));
Xstar := CreateShimuraQuot(D, N, Set(Divisors(D*N)));
Xstar`g := GenusShimuraCurveQuotient(D, N, Xstar`W); Xstar`CurveID := 0;
curves := GetQuotientsAndGenera([Xstar]);
_ := exists(star){c : c in curves | IsStarCurve(c)};
Ldata := ShimuraCurveLattice(D, N);
emit(Sprintf("det Q = %o = %o, M = %o", Determinant(ChangeRing(Ldata`Q, Integers())),
             Factorization(Determinant(ChangeRing(Ldata`Q, Integers()))),
             IsOdd(D*N) select 4*D*N else 2*D*N));
fs := BorcherdsForms(star, curves : Prec := 100);
ks := Sort([k : k in Keys(fs)]);
emit(Sprintf("keys = %o", ks));
for k in ks do
    emit(Sprintf("--------- form[%o] ---------", k));
    dump_form(qExpansionAtoo(fs[k], 1), qExpansionAt0(fs[k], 1), Ldata, D, N);
end for;
emit("DONE");
exit;
