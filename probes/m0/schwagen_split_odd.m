// Gap-1 validation across bases, on the FULL BorcherdsForms list (not just fs[-1]).
// Compares the OLD local factor (LocalWhittakerAtOne) against the NEW Schwagenscheidt Thm 1.4
// density N_{gamma,n}(p^{w_p})/p^{2w_p} with Eisenstein Gram = +Q_code (settled empirically on 15_2).
// PASS criteria: 15_2 -> 4 on the hauptmodul and integers on every form; 21_2 -> integers (old: -20/3).

AttachSpec("ShimuraQuotients.spec");
fh := Open("split_odd_out.txt", "w");
emit := procedure(s) fprintf fh, "%o\n", s; Flush(fh); end procedure;

rep_count := function(G, gam, n, p, w)
    pw := p^w;
    if w le 0 then return 1, pw; end if;
    Gr := ChangeRing(G, Rationals());
    Gg := gam*Gr;
    cst := (gam*Gr, gam)/2 + n;
    assert IsIntegral(cst);
    for i in [1..3] do assert IsIntegral(Gg[i]); end for;
    h11 := Integers()!(G[1][1]/2); h22 := Integers()!(G[2][2]/2); h33 := Integers()!(G[3][3]/2);
    g12 := Integers()!G[1][2]; g13 := Integers()!G[1][3]; g23 := Integers()!G[2][3];
    b1 := Integers()!Gg[1]; b2 := Integers()!Gg[2]; b3 := Integers()!Gg[3];
    c0 := Integers()!cst;
    Am := h33 mod pw;
    hist := [];
    for b in [0..pw-1] do
        row := [0 : v in [1..pw]];
        for t in [0..pw-1] do
            v := (Am*t*t + b*t) mod pw;
            row[v+1] +:= 1;
        end for;
        Append(~hist, row);
    end for;
    N := 0;
    for x1 in [0..pw-1] do
        for x2 in [0..pw-1] do
            B := (g13*x1 + g23*x2 - b3) mod pw;
            C := (h11*x1*x1 + h22*x2*x2 + g12*x1*x2 - b1*x1 - b2*x2 + c0) mod pw;
            N +:= hist[B+1][((-C) mod pw)+1];
        end for;
    end for;
    return N, pw;
end function;

// mode 0 = OLD (LocalWhittakerAtOne, bucket +j)
// mode 1 = density, Eisenstein Gram +Q_code  => Thm 1.4 forces bucket -j
// mode 2 = density, Eisenstein Gram -Q_code  => Thm 1.4 forces bucket +j
// (modes 1,2 are the only two combinations consistent with n in Z - Q(gamma).)
multiplier := function(foo, f0, Ldata, D, N, mode)
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
        // N_gamma = order of gamma in L'/L = smallest k>0 with k*gamma in L = Z^3
        Ngam := LCM([Denominator(x) : x in Eltseq(w_eta)]);
        for p in Sc do
            if mode ne 0 then
                Guse := (mode eq 1) select Qint else negQ;
                w := 1 + 2*Valuation(Rationals()!(2*Ngam*r), p);
                cnt, pw := rep_count(Guse, Vector(Rationals(), Eltseq(w_eta)), r, p, w);
                loc := (Rationals()!cnt)/(Rationals()!pw)^2;
            else
                loc := LocalWhittakerAtOne(r, p, Vector(Rationals(), Eltseq(w_eta)), Lfull, negQ);
            end if;
            g *:= loc;
            if g eq 0 then break; end if;
        end for;
        if g eq 0 then return Rationals()!0; end if;
        en := &*[Rationals() | 1 - Evaluate(chi, p)/p : p in Sc];
        ed := &*[Rationals() | 1 - 1/(Rationals()!p)^2 : p in Sc];
        return c * cond_half * (Rationals()!h/wr) * (en/ed) * g;
    end function;

    Too := Rationals()!0; T0 := Rationals()!0;
    T := Rationals()!0;
    for mm in [1..-Valuation(foo)] do
        c := Coefficient(foo, -mm);
        if c ne 0 then Too +:= contrib(i0, Rationals()!mm, Rationals()!c); end if;
    end for;
    for j in [1..-Valuation(f0)] do
        c := Coefficient(f0, -j);
        if c eq 0 then continue; end if;
        r := (Rationals()!j)/M;
        // Thm 1.4 requires n in Z - Q(gamma); with Eisenstein Gram +Q_code that means
        // r = -nm_eta mod 1, i.e. the bucket at -j (not +j). Densities are insensitive to
        // gamma -> -gamma (send r -> -r in the count), so only the pairing sign matters.
        bucket := (mode eq 1) select ((-j) mod M) else (j mod M);
        for eta in mod_M_to_vecs[bucket] do T0 +:= contrib(eta, r, Rationals()!c); end for;
    end for;
    return -96 * (Too+T0) / (D*N), -96*Too/(D*N), -96*T0/(D*N);
end function;

emit("START");
for base in [<10,11>] do
    D, N := Explode(base);
    emit(Sprintf("=========== X0^%o(%o) ===========", D, N));
    try
        Xstar := CreateShimuraQuot(D, N, Set(Divisors(D*N)));
        Xstar`g := GenusShimuraCurveQuotient(D, N, Xstar`W); Xstar`CurveID := 0;
        curves := GetQuotientsAndGenera([Xstar]);
        _ := exists(star){c : c in curves | IsStarCurve(c)};
        Ldata := ShimuraCurveLattice(D, N);
        fs := BorcherdsForms(star, curves : Prec := 100);
        ks := Sort([k : k in Keys(fs)]);
        emit(Sprintf("  #forms = %o, keys = %o", #ks, ks));
        for k in ks do
            eta := fs[k];
            foo := qExpansionAtoo(eta, 1);
            f0  := qExpansionAt0(eta, 1);
            m0, o0, z0 := multiplier(foo, f0, Ldata, D, N, 0);
            m1, o1, z1 := multiplier(foo, f0, Ldata, D, N, 1);
            emit(Sprintf("  form[%-3o]: OLD tot=%-8o (oo=%-8o zero=%-8o) | dens tot=%-8o (oo=%-8o zero=%-8o)",
                         k, m0, o0, z0, m1, o1, z1));
        end for;
    catch e
        emit(Sprintf("  ERROR building base: %o", e`Object));
    end try;
end for;
emit("DONE");
