// ctThetaAt: ctTheta with an ARBITRARY tracked coset estc, for N > 1 where
// the tracked coset is a NONZERO isotropic element.  Same closed form; only
// the target of c*y = tgt, the T-coset value and an isotropy assertion move.
ctThetaAt := function(Gram, words, CC, QSIGN, estc)
    ii := CC.1; pi := Pi(CC);
    ee := func< z | Exp(2*pi*ii*z) >;
    frac := func< r | r - Floor(r) >;

    dn := Denominator(Gram^(-1));
    Ldual := RSpaceWithBasis(ChangeRing(dn*Gram^(-1), Integers()));
    Lp := RSpaceWithBasis(ScalarMatrix(3, dn));
    G, toG := Ldual/Lp;

    mods := Moduli(G);
    keep := [ r : r in [1..#mods] | mods[r] gt 1 ];
    ms := [ mods[r] : r in keep ];
    k := #ms;
    n := &*ms;

    gens := [ G.(keep[r]) : r in [1..k] ];
    wg := [ ChangeRing(g@@toG, Rationals()) : g in gens ];
    QG := [ frac(QSIGN*(wg[r]*Gram, wg[r])/(2*dn^2)) : r in [1..k] ];
    BG := [ [ frac(QSIGN*(wg[r]*Gram, wg[s])/dn^2) : s in [1..k] ] : r in [1..k] ];

    zero := [ Integers() | 0 : r in [1..k] ];
    mulc := func< t, c | [ (t*c[r]) mod ms[r] : r in [1..k] ] >;
    iszc := func< c | forall{ r : r in [1..k] | c[r] mod ms[r] eq 0 } >;
    Qc := function(c)
        s := &+[ Rationals() | c[r]^2*QG[r] : r in [1..k] ];
        for r in [1..k] do for s2 in [r+1..k] do s +:= c[r]*c[s2]*BG[r][s2]; end for; end for;
        return frac(s);
    end function;
    Bc := function(c, e)
        s := Rationals()!0;
        for r in [1..k] do for s2 in [1..k] do s +:= c[r]*e[s2]*BG[r][s2]; end for; end for;
        return frac(s);
    end function;

    levD := Lcm([ Integers() | Denominator(QG[r]) : r in [1..k] ]
            cat [ Integers() | Denominator(BG[r][s2]) : r, s2 in [1..k] ]);

    plist := PrimeDivisors(n);
    JP := AssociativeArray(); LEVP := AssociativeArray();
    for p in plist do
        cur := [ zero ];
        for r in [1..k] do
            ar := Valuation(ms[r], p); step := ms[r] div p^ar;
            nxt := [ ];
            for c in cur do
                for t in [0..p^ar - 1] do
                    cc := c; cc[r] := (t*step) mod ms[r]; Append(~nxt, cc);
                end for;
            end for;
            cur := nxt;
        end for;
        JP[p] := cur;
        ex := &*[ Integers() | p^Valuation(ms[r], p) : r in [1..k] ];
        LEVP[p] := LCM([ Integers() | ex ] cat [ Denominator(Qc(j)) : j in cur ]);
    end for;

    // Jordan basis of the 2-part and the canonical x_c (Stromberg Def 2.15)
    xc := zero; r2 := 0;
    if 2 in plist then
        J2 := JP[2];
        tor2 := [ j : j in J2 | iszc(mulc(2, j)) ];
        tor2nz := [ j : j in tor2 | not iszc(j) ];
        r2 := Valuation(#tor2, 2);
        error if #tor2nz ne 2^r2 - 1, "2-torsion count mismatch";
        error if r2 notin {1, 3}, "unexpected 2-rank in a Gross lattice", r2;
        basis2 := [ ]; found := false;
        if r2 eq 1 then
            basis2 := [ tor2nz[1] ]; found := true;
        end if;
        for i1 in (found select [ ] else tor2nz) do
            for i2 in tor2nz do
                if i2 eq i1 or Bc(i1, i2) ne 0 then continue; end if;
                for i3 in tor2nz do
                    if i3 eq i1 or i3 eq i2 then continue; end if;
                    if i3 eq [ (i1[r] + i2[r]) mod ms[r] : r in [1..k] ] then continue; end if;
                    if Bc(i1, i3) ne 0 or Bc(i2, i3) ne 0 then continue; end if;
                    basis2 := [i1, i2, i3]; found := true; break;
                end for;
                if found then break; end if;
            end for;
            if found then break; end if;
        end for;
        error if not found, "no Jordan basis for J_2";
        for b in basis2 do xc := [ (xc[r] + b[r]) mod ms[r] : r in [1..k] ]; end for;
    end if;

    GTAB := AssociativeArray();
    for p in plist do
        J := JP[p];
        QJp := [ Qc(j) : j in J ]; BXCp := [ Bc(xc, j) : j in J ];
        row := [ ];
        for rr in [0..LEVP[p] - 1] do
            S0c := CC!0; S1c := CC!0; nk := 0;
            for jx in [1..#J] do
                t := rr*QJp[jx];
                S0c +:= ee(frac(t)); S1c +:= ee(frac(t + BXCp[jx]));
                if iszc(mulc(rr, J[jx])) then nk +:= 1; end if;
            end for;
            nrm := Sqrt(CC!(#J) * nk);
            Append(~row, <S0c/nrm, S1c/nrm>);
        end for;
        GTAB[p] := row;
    end for;
    gsumJ := func< r, useXc, p |
        useXc select GTAB[p][(r mod LEVP[p]) + 1][2] else GTAB[p][(r mod LEVP[p]) + 1][1] >;

    R30 := RealField(30);
    sig8 := func< p | Round(8*Arg(gsumJ(1, false, p))/(2*Pi(R30))) mod 8 >;
    signD := (&+[ Integers() | sig8(p) : p in plist ]) mod 8;
    sign2 := 2 in plist select sig8(2) else 0;

    hilb := func< x, y | (x lt 0 and y lt 0) select -1 else 1 >;
    oddpart := function(c) c2 := AbsoluteValue(c); while IsEven(c2) do c2 := c2 div 2; end while; return c2; end function;
    sfun := func< g | g[2][1] ne 0 select g[2][1] else g[2][2] >;
    Smat := Matrix(Integers(), 2, 2, [0, -1, 1, 0]);
    Tmat := func< kk | Matrix(Integers(), 2, 2, [1, kk, 0, 1]) >;
    kubota := function(w)
        ls := [ ];
        for t in w do Append(~ls, t[1] eq "S" select Smat^(-1) else Tmat(-t[2])); end for;
        acc := ls[1]; sign := 1;
        for i in [2..#ls] do
            h := ls[i]; prod := h*acc;
            sign *:= (sfun(h)*sfun(prod) lt 0 and sfun(acc)*sfun(prod) lt 0) select -1 else 1;
            acc := prod;
        end for;
        return sign, acc;
    end function;

    out := [ CC | 0 : w in words ];
    YTab := AssociativeArray();
    for wi->w in words do
        if #w eq 0 or VVWordMatrix(w)[2][1] eq 0 then
            out[wi] := iszc(estc) select CC!1 else CC!0;   // T^k: e_0 read at estc
            continue;
        end if;
        gmat := VVWordMatrix(w); gi := gmat^(-1);
        a := gi[1][1]; c := gi[2][1];
        if a eq 0 then
            // The S-coset constant is e(sig/8)/sqrt|A| where sig is the
            // SIGNATURE of the discriminant form mod 8, not the universal
            // e(1/8).  eis32s hardcodes e(1/8) and is right there only because
            // every Shimura-curve lattice in this campaign reports
            // sign(Dbar) = 1 (it uses the NEGATED signature-(1,2) form, so
            // sig = +1).  With QSIGN = -1 these Gross lattices are NEGATIVE
            // definite, sig = -3 = 5 mod 8, and e(5/8) = -e(1/8): a sign flip.
            // Diagnosed at 330_1, where hardcoding e(1/8) left exactly ONE
            // coset of 3456 wrong, by |4 v1 - 6 v2| = 0.004288 predicted
            // against 0.0042855 measured.  signD is what the Gauss sums say,
            // so use it.
            error if #[ t : t in w | t[1] eq "S" ] ne 1,
                "a = 0 word has more than one S", wi;
            error if Qc(estc) ne 0, "tracked coset not isotropic";
            out[wi] := ee(CC!signD/8) / Sqrt(CC!n);
            continue;
        end if;
        useXc := (Valuation(c, 2) eq 1);
        cm := ((c mod levD) + levD) mod levD;
        key := <cm, useXc>;
        if not IsDefined(YTab, key) then
            tgt := useXc select [ (estc[r] - xc[r]) mod ms[r] : r in [1..k] ] else estc;
            okY := true; y := [ Integers() | ];
            for r in [1..k] do
                m := ms[r]; g0 := GCD(cm, m);
                if tgt[r] mod g0 ne 0 then okY := false; break; end if;
                mm := m div g0;
                if mm eq 1 then Append(~y, 0);
                else Append(~y, (InverseMod((cm div g0) mod mm, mm) * ((tgt[r] div g0) mod mm)) mod mm);
                end if;
            end for;
            dcn := &*[ Integers() | GCD(cm, ms[r]) : r in [1..k] ];
            YTab[key] := <okY, y, dcn>;
        end if;
        okY := YTab[key][1]; y := YTab[key][2]; dcn := YTab[key][3];
        if not okY then continue; end if;
        xiJ := CC!1;
        for p in plist do
            if c mod p ne 0 then xiJ *:= gsumJ(c, false, p);
            else xiJ *:= CC!KroneckerSymbol(-a, #JP[p]) * gsumJ(-a*c, useXc, p);
            end if;
        end for;
        xi0 := CC!KroneckerSymbol(-a, AbsoluteValue(c)) * hilb(-a, c);
        xi2 := IsOdd(c) select CC!1 else ee(CC!(-(a + 1)*(oddpart(c) - 1 + sign2))/8);
        xidef := ee(CC!(-signD)/4) * xi0 * xi2 * xiJ;
        ksign, acc := kubota(w);
        error if acc ne gi and acc ne -gi, "word matrix mismatch", wi;
        u := ksign;
        if #w eq 1 and w[1][1] eq "S" then u := -u; end if;
        out[wi] := u * xidef * Sqrt(CC!dcn/CC!n) *
            ee(frac(a*c*Qc(y) + (useXc select Bc(xc, y) else 0)));
    end for;
    return out, n, ms;
end function;
