// N > 1: LOCALISE the class-dependent residue.
//
// The N=1 theorem says ct factors over primes.  At N>1 the Shimura lattice is
// Eichler (hyperbolic) at p | N and the tracked coset mu is NONZERO there, so
//     rho_mu = [ prod_{p|D} c_p^ram ] * h_N(mu) * (2-adic extras)
// while the Gross combination pred0 = 1/2 two - 1/2 v3 carries c_N^Eich or
// c_N^ram in that slot.  The measured residue was lambda = 2N/(N+1) on most
// cusp classes and 1 on the classes {q,4,8,qr,4r,8r}.  This script separates
// the pieces that can carry it:
//     kappa   = rho_mu / rho_0     (the OFF-DIAGONAL / DIAGONAL ratio: this is
//                                   the only mu-dependent factor, and the
//                                   closed form says |kappa| in {0,1})
//     base    = rho_0 / pred0      (an N=1-style comparison, no mu at all)
// If |kappa| is 0 or 1 everywhere then the magnitude 4/3 vs 1 CANNOT come from
// the tracked coset, and the residue lives entirely in base.
AttachSpec("ShimuraQuotients.spec");
load "vvdata/weyl-campaign/fastcosets.m";
load "vvdata/weyl-campaign/etared.m";
load "vvdata/weyl-campaign/ctlat.m";
load "vvdata/weyl-campaign/ctlatat.m";
PREC := 120; CC := ComplexField(PREC); QSIGN := -1;
frac := func< r | r - Floor(r) >;
mass := function(Dps, Rs)
    n := &*[ Integers() | p-1 : p in Dps ];
    m := Rationals()!n / (48*2^(#Dps-1));
    for x in Rs do m *:= Rationals()!(x+1)/2; end for;
    return m;
end function;
R8 := ComplexField(8);
RR := RealField(8);
tol := RealField(PREC)!10^(-80);

//        D,  N, D', R, Rs
cases := [ <22,3,11,3,6>, <15,2,5,2,6>, <10,7,5,7,14> ];

for cs in cases do
    D := cs[1]; N := cs[2]; Dp := cs[3]; Rl := cs[4]; Rls := cs[5];
    DN := D*N; lev := 4*DN;
    words := [ VVSTWord(g) : g in fastCosetReps(lev) ];
    nw := #words;
    Ld := ShimuraCurveLattice(D, N);
    Gm := ChangeRing(Ld`Q, Rationals()); dn := Ld`denom;
    Gg := Ld`disc_grp; elts := [ g : g in Gg ];
    mods := Moduli(Gg); keep := [r : r in [1..#mods] | mods[r] gt 1];
    ms_ := [mods[r] : r in keep]; k := #ms_;
    Qb := [ frac(-(ChangeRing(e@@Ld`to_disc, Rationals())*Gm,
                   ChangeRing(e@@Ld`to_disc, Rationals()))/(2*dn^2)) : e in elts ];
    iso := [ i : i in [1..#elts] | Qb[i] eq 0 ];
    estc := [ (Eltseq(elts[iso[2]])[keep[r]]) mod ms_[r] : r in [1..k] ];
    ord := LCM([ Integers() | ms_[r] div GCD(ms_[r], estc[r]) : r in [1..k] ]);
    printf "\n=== %o_%o   DN=%o  lev=%o  |A|=%o  moduli=%o  #iso=%o\n",
        D, N, DN, lev, &*ms_, ms_, #iso;
    printf "    tracked coset estc = %o,  order %o\n", estc, ord;

    rho  := ctThetaAt(Gm, words, CC, QSIGN, estc);
    rho0 := ctTheta(Gm, words, CC, QSIGN);
    m1 := mass([Dp],[Rl]); msx := mass([Dp],[Rls]);
    v1 := ctTheta(grossGram(Dp, Rl), words, CC, QSIGN);
    vs := ctTheta(grossGram(Dp, Rls), words, CC, QSIGN);
    v3 := ctTheta(grossGram(DN, 1), words, CC, QSIGN);
    two := [ CC | (CC!msx*vs[i] - CC!m1*v1[i])/CC!(msx-m1) : i in [1..nw] ];
    pr  := [ CC | (CC!(1/2))*two[i] - (CC!(1/2))*v3[i] : i in [1..nw] ];

    cls := AssociativeArray();
    for i in [1..nw] do
        g := VVWordMatrix(words[i]); c := g[2][1];
        gg := GCD(c mod lev, lev);
        if not IsDefined(cls, gg) then cls[gg] := [Integers()|]; end if;
        Append(~cls[gg], i);
    end for;

    setof := function(idx, num, den)
        vals := {@ @};
        for i in idx do
            if Abs(den[i]) gt tol then Include(~vals, R8!(num[i]/den[i]));
            elif Abs(num[i]) gt tol then Include(~vals, R8!(10^9));
            end if;
        end for;
        return [ vals[j] : j in [1..Minimum(3,#vals)] ], #vals;
    end function;

    printf "  %-6o %-5o %-4o %-4o %-4o | %-26o | %o\n",
        "class", "#", "nz0", "nzM", "nzP", "kappa = rho_mu/rho_0", "base = rho_0/pred0";
    for gg in Sort([x : x in Keys(cls)]) do
        idx := cls[gg];
        nz0 := #[ i : i in idx | Abs(rho0[i]) gt tol ];
        nzM := #[ i : i in idx | Abs(rho[i])  gt tol ];
        nzP := #[ i : i in idx | Abs(pr[i])   gt tol ];
        kv, nk := setof(idx, rho, rho0);
        bv, nb := setof(idx, rho0, pr);
        printf "  %-6o %-5o %-4o %-4o %-4o | %-26o | %o\n",
            gg, #idx, nz0, nzM, nzP,
            Sprintf("%o%o", kv, nk gt 3 select Sprintf("+%o", nk-3) else ""),
            Sprintf("%o%o", bv, nb gt 3 select Sprintf("+%o", nb-3) else "");
    end for;

    // the full residue, for the record
    printf "  --- residue rho_mu/pred0 by class ---\n";
    for gg in Sort([x : x in Keys(cls)]) do
        rv, nr := setof(cls[gg], rho, pr);
        printf "  %-6o %o%o\n", gg, rv, nr gt 3 select Sprintf("+%o", nr-3) else "";
    end for;
end for;
printf "DONE\n";
quit;
