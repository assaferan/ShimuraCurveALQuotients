// N > 1: test the conjectured three-term identity, BOTH SIDES at the FORM level 4DN.
//    rho_est =?= 1/2 [Theta(D',Rs) - Theta(D',R)]/(mass_Rs - mass_R) - 1/2 thetabar(DN,1)
// rho_est via ctThetaAt on the Shimura lattice at the tracked (nonzero) coset;
// the theta side via ctTheta (its (0,0) entries).  Same coset list throughout.
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

cases := [ <15, 2, 5, 2, 6>, <22, 3, 11, 3, 6> ];
for c in cases do
    D := c[1]; N := c[2]; Dp := c[3]; Rl := c[4]; Rls := c[5];
    DN := D*N; lev := 4*DN;
    printf "\n===== X_0^%o(%o), DN = %o, form level %o =====\n", D, N, DN, lev;
    words := [ VVSTWord(g) : g in fastCosetReps(lev) ];
    nw := #words;

    Ld := ShimuraCurveLattice(D, N);
    Gm := ChangeRing(Ld`Q, Rationals()); dn := Ld`denom;
    Gg := Ld`disc_grp; elts := [ g : g in Gg ];
    mods := Moduli(Gg); keep := [r : r in [1..#mods] | mods[r] gt 1];
    ms_ := [mods[r] : r in keep]; k := #ms_;
    coord := func< e | [ (Eltseq(e)[keep[r]]) mod ms_[r] : r in [1..k] ] >;
    Qb := [ frac(-(ChangeRing(e@@Ld`to_disc, Rationals())*Gm,
                   ChangeRing(e@@Ld`to_disc, Rationals()))/(2*dn^2)) : e in elts ];
    iso := [ i : i in [1..#elts] | Qb[i] eq 0 ];
    estc := coord(elts[iso[2]]);
    printf "  #isotropic = %o, tracked coset coords %o\n", #iso, estc;

    rho := ctThetaAt(Gm, words, CC, QSIGN, estc);
    rn := Maximum([ Abs(x) : x in rho ]);
    printf "  |rho_est|_max = %o\n", RealField(8)!rn;

    m1 := mass([Dp], [Rl]); msx := mass([Dp], [Rls]);
    v1 := ctTheta(grossGram(Dp, Rl), words, CC, QSIGN);
    vs := ctTheta(grossGram(Dp, Rls), words, CC, QSIGN);
    v3 := ctTheta(grossGram(DN, 1), words, CC, QSIGN);
    pred := [ CC | (CC!(1/2))*((CC!msx*vs[i] - CC!m1*v1[i])/CC!(msx-m1))
                   - (CC!(1/2))*v3[i] : i in [1..nw] ];
    dev := Maximum([ Abs(pred[i] - rho[i]) : i in [1..nw] ]);
    printf "  conjecture 1/2[Theta(%o,%o)-Theta(%o,%o)]/(%o) - 1/2 thetabar(%o,1)\n",
        Dp, Rls, Dp, Rl, msx-m1, DN;
    printf "  max|conjecture - rho_est| = %o  (rel %o)  %o\n",
        RealField(8)!dev, RealField(8)!(dev/rn),
        dev lt 10^(-90)*rn select "*** HOLDS ***" else "does not hold as stated";
    if dev ge 10^(-90)*rn then
        // is the discrepancy a scalar multiple of one of the pieces?
        printf "    (rho_est/pred at a few cosets: %o)\n",
            [ ComplexField(8)!(Abs(pred[i]) gt 10^(-60) select rho[i]/pred[i] else CC!0)
              : i in [2..5] ];
    end if;
end for;
printf "DONE\n";
quit;
