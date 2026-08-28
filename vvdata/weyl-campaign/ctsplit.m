// ctSplitLocal: the PER-PRIME factorisation of the constant-term entry
//
//     ct_L(gamma)  =  glob(gamma) * prod_p  loc_p^L(gamma)
//
// with glob = u * xi0 depending ONLY on gamma (the Kubota sign and the
// Kronecker/Hilbert symbol), so it is the SAME for every lattice.  Everything
// lattice-dependent is local:
//
//     loc_p = e(-sig8(p)/4) * [ gauss sum at p ] * sqrt( dcn_p / n_p )
//
// and at p = 2 also the xi2 factor, which depends on the lattice only through
// sign2 = sig8(2).  This is the tool for comparing two lattices prime by
// prime: ct_L / ct_G = prod_p loc_p^L / loc_p^G, no global phase left over.
//
// RESTRICTION.  Written for the ZERO tracked coset and for cosets with
// ord_2(c) != 1.  Those are exactly the two places the closed form needs the
// Stromberg x_c correction and the y-phase, and both are trivial otherwise
// (y = 0 solves c*y = 0).  Words with c = 0 or a = 0 are returned as invalid.
ctSplitLocal := function(Gram, words, CC, QSIGN)
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

    levD := Lcm([ Integers() | Denominator(QG[r]) : r in [1..k] ]
            cat [ Integers() | Denominator(BG[r][s2]) : r, s2 in [1..k] ]);

    plist := PrimeDivisors(n);
    JP := AssociativeArray(); LEVP := AssociativeArray(); NP := AssociativeArray();
    MSP := AssociativeArray();
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
        NP[p] := #cur;
        MSP[p] := [ p^Valuation(ms[r], p) : r in [1..k] ];
        ex := &*[ Integers() | p^Valuation(ms[r], p) : r in [1..k] ];
        LEVP[p] := LCM([ Integers() | ex ] cat [ Denominator(Qc(j)) : j in cur ]);
    end for;

    GTAB := AssociativeArray();
    for p in plist do
        J := JP[p];
        QJp := [ Qc(j) : j in J ];
        row := [ ];
        for rr in [0..LEVP[p] - 1] do
            S0c := CC!0; nk := 0;
            for jx in [1..#J] do
                S0c +:= ee(frac(rr*QJp[jx]));
                if iszc(mulc(rr, J[jx])) then nk +:= 1; end if;
            end for;
            Append(~row, S0c/Sqrt(CC!(#J) * nk));
        end for;
        GTAB[p] := row;
    end for;
    gsum := func< r, p | GTAB[p][(r mod LEVP[p]) + 1] >;

    R30 := RealField(30);
    sig8 := func< p | Round(8*Arg(gsum(1, p))/(2*Pi(R30))) mod 8 >;
    SIG := AssociativeArray();
    for p in plist do SIG[p] := sig8(p); end for;
    sign2 := 2 in plist select SIG[2] else 0;

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

    OK := [ Booleans() | ]; GLOB := [ CC | ]; LOC := [* *];
    for wi->w in words do
        ok := true;
        if #w eq 0 or VVWordMatrix(w)[2][1] eq 0 then ok := false; end if;
        if ok then
            gmat := VVWordMatrix(w); gi := gmat^(-1);
            a := gi[1][1]; c := gi[2][1];
            if a eq 0 or Valuation(c, 2) eq 1 then ok := false; end if;
        end if;
        if not ok then
            Append(~OK, false); Append(~GLOB, CC!0); Append(~LOC, AssociativeArray());
            continue;
        end if;
        gi := VVWordMatrix(w)^(-1);
        a := gi[1][1]; c := gi[2][1];
        cm := ((c mod levD) + levD) mod levD;
        lp := AssociativeArray();
        for p in plist do
            dcnp := &*[ Integers() | GCD(cm, MSP[p][r]) : r in [1..k] ];
            if c mod p ne 0 then g := gsum(c, p);
            else g := CC!KroneckerSymbol(-a, NP[p]) * gsum(-a*c, p);
            end if;
            f := ee(CC!(-SIG[p])/4) * g * Sqrt(CC!dcnp/CC!NP[p]);
            if p eq 2 then
                f *:= IsOdd(c) select CC!1
                      else ee(CC!(-(a + 1)*(oddpart(c) - 1 + sign2))/8);
            end if;
            lp[p] := f;
        end for;
        ksign, acc := kubota(w);
        error if acc ne gi and acc ne -gi, "word matrix mismatch", wi;
        u := ksign;
        if #w eq 1 and w[1][1] eq "S" then u := -u; end if;
        Append(~OK, true);
        Append(~GLOB, CC!u * CC!KroneckerSymbol(-a, AbsoluteValue(c)) * hilb(-a, c));
        Append(~LOC, lp);
    end for;
    return OK, GLOB, LOC, plist, ms, SIG;
end function;
