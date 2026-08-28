// N > 1: the three-term identity, with the scalar pinned.
//   rho_est  =?=  lambda(N) * { 1/2 [Theta(D',Rs)-Theta(D',R)]/(mass_Rs-mass_R)
//                               - 1/2 thetabar(DN,1) },   lambda(N) = 2N/(N+1).
// 2N/(N+1) = N / [(N+1)/2], and (N+1)/2 is the LOCAL MASS FACTOR at the Eichler
// prime N -- which is why this is the natural normalisation to try.
// 10_7 (N=7, lambda = 7/4) is the discriminating case.
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

//        D,  N, D', R, Rs
cases := [ <15,2,5,2,6>, <21,2,7,2,6>, <22,3,11,3,6>, <33,2,11,2,6>,
           <35,2,7,2,10>, <55,2,11,2,10>, <10,7,5,7,14> ];
printf "%-8o %-4o %-10o %-16o %o\n", "base", "N", "2N/(N+1)", "measured", "verdict";
for c in cases do
    D := c[1]; N := c[2]; Dp := c[3]; Rl := c[4]; Rls := c[5];
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
    rho := ctThetaAt(Gm, words, CC, QSIGN, estc);
    m1 := mass([Dp], [Rl]); msx := mass([Dp], [Rls]);
    v1 := ctTheta(grossGram(Dp, Rl), words, CC, QSIGN);
    vs := ctTheta(grossGram(Dp, Rls), words, CC, QSIGN);
    v3 := ctTheta(grossGram(DN, 1), words, CC, QSIGN);
    lam := CC!(2*N)/CC!(N+1);
    pred := [ CC | lam*((CC!(1/2))*((CC!msx*vs[i] - CC!m1*v1[i])/CC!(msx-m1))
                        - (CC!(1/2))*v3[i]) : i in [1..nw] ];
    rn := Maximum([ Abs(x) : x in rho ]);
    dev := Maximum([ Abs(pred[i] - rho[i]) : i in [1..nw] ]);
    // also report the raw measured ratio, for the record
    rat := CC!0;
    for i in [1..nw] do
        p0 := (CC!(1/2))*((CC!msx*vs[i] - CC!m1*v1[i])/CC!(msx-m1)) - (CC!(1/2))*v3[i];
        if Abs(p0) gt 10^(-40) then rat := rho[i]/p0; break; end if;
    end for;
    printf "%o_%-6o N=%o lambda=%o measured=%o %o\n", D, N, N,
        RationalField()!(2*N/(N+1)), ComplexField(10)!rat,
        dev lt 10^(-90)*rn select "HOLDS" else Sprintf("MISMATCH rel %o", RealField(6)!(dev/rn));
    if dev ge 10^(-90)*rn then
        // WHICH cosets disagree, and how do they sit by cusp class?
        okc := AssociativeArray(); badc := AssociativeArray();
        for i in [1..nw] do
            g := VVWordMatrix(words[i]); cls := GCD(g[2][1] mod lev, lev);
            good := Abs(pred[i] - rho[i]) lt 10^(-90)*rn;
            if good then okc[cls] := (IsDefined(okc,cls) select okc[cls] else 0) + 1;
            else badc[cls] := (IsDefined(badc,cls) select badc[cls] else 0) + 1; end if;
        end for;
        ks := Sort([ c : c in Keys(badc) ]);
        printf "    disagree on %o of %o cosets; by cusp class gcd(c,%o): %o\n",
            &+[ badc[c] : c in ks ], nw, lev, [ <c, badc[c]> : c in ks ];
        printf "    agree by class: %o\n",
            [ <c, okc[c]> : c in Sort([ x : x in Keys(okc) ]) ];
        // on the failing cosets, is rho/pred a DIFFERENT constant?
        vals := {@ @};
        for i in [1..nw] do
            if Abs(pred[i] - rho[i]) ge 10^(-90)*rn and Abs(pred[i]) gt 10^(-40) then
                Include(~vals, ComplexField(8)!(rho[i]/pred[i]));
            end if;
        end for;
        printf "    rho/pred on the failing cosets: %o distinct, first few %o\n",
            #vals, [ vals[j] : j in [1..Minimum(4,#vals)] ];
    end if;
end for;
printf "DONE\n";
quit;
