// N > 1: test the conjectured three-term identity
//    rho_est = 1/2 [Theta(D',Rs) - Theta(D',R)]/(mass_Rs - mass_R) - 1/2 thetabar(DN,1)
// read off the banked fits (15_2, 22_3), where the two-term part is exactly
// HALF the N=1 weights and the third term is always -1/2 thetabar(DN,1).
//
// rho_est is the Weil entry at a NONZERO tracked coset, so it is NOT the (0,0)
// entry the N=1 theorem uses.  Here it is taken from eis32s (which handles est
// properly); the theta side is (0,0) entries via ctTheta.
AttachSpec("ShimuraQuotients.spec");
load "vvdata/weyl-campaign/fastcosets.m";
load "vvdata/weyl-campaign/etared.m";
load "vvdata/weyl-campaign/ctlat.m";
PREC := 120; CC := ComplexField(PREC); QSIGN := -1;

mass := function(Dps, Rs)
    n := &*[ Integers() | p-1 : p in Dps ];
    m := Rationals()!n / (48*2^(#Dps-1));
    for x in Rs do m *:= Rationals()!(x+1)/2; end for;
    return m;
end function;

// D, N, and the banked support (D', R, Rs) plus the third slot (DN,1)
cases := [ <15, 2, 5, [2], [2,3]>, <22, 3, 11, [3], [3,2]> ];
for c in cases do
    D := c[1]; N := c[2]; Dp := c[3]; Rl := c[4]; Rls := c[5];
    DN := D*N; lev := 4*DN;
    printf "\n===== X_0^%o(%o), DN = %o, level %o =====\n", D, N, DN, lev;
    words := [ VVSTWord(g) : g in fastCosetReps(lev) ];
    nw := #words;

    // rho at the TRACKED (nonzero) coset, straight from the Shimura lattice
    Ld := ShimuraCurveLattice(D, N);
    // reproduce eis32s's est: second isotropic coset in the library's order
    G := Ld`disc_grp; elts := [ g : g in G ];
    Qr := ChangeRing(Ld`Q, Rationals()); dn := Ld`denom;
    frac := func< r | r - Floor(r) >;
    Qb := [ frac(-(ChangeRing(e@@Ld`to_disc, Rationals())*Qr,
                   ChangeRing(e@@Ld`to_disc, Rationals()))/(2*dn^2)) : e in elts ];
    iso := [ i : i in [1..#elts] | Qb[i] eq 0 ];
    printf "  #isotropic = %o (2N-1 = %o); tracked coset index %o\n", #iso, 2*N-1, iso[2];

    // theta pieces
    m1 := mass([Dp], Rl); ms := mass([Dp], Rls);
    v1 := ctTheta(grossGram(Dp, &*Rl), words, CC, QSIGN);
    vs := ctTheta(grossGram(Dp, &*Rls), words, CC, QSIGN);
    v3 := ctTheta(grossGram(DN, 1), words, CC, QSIGN);
    printf "  mass(%o,%o) = %o, mass(%o,%o) = %o, difference %o\n",
        Dp, &*Rl, m1, Dp, &*Rls, ms, ms-m1;
    pred := [ CC | (CC!(1/2))*((CC!ms*vs[i] - CC!m1*v1[i])/CC!(ms-m1)) - (CC!(1/2))*v3[i]
              : i in [1..nw] ];
    printf "  conjecture: 1/2 [Theta(%o,%o)-Theta(%o,%o)]/(mass diff) - 1/2 thetabar(%o,1)\n",
        Dp, &*Rls, Dp, &*Rl, DN;
    printf "  predicted a(0) at the identity coset = %o   (weights sum to 0, so expect 0)\n",
        ComplexField(12)!pred[1];
    // what does rho look like there?
    printf "  NOTE: rho_est must come from eis32s (nonzero coset); ctTheta gives (0,0) only.\n";
    printf "  |predicted|_max = %o\n", RealField(8)!Maximum([Abs(x) : x in pred]);
end for;
printf "DONE\n";
quit;
