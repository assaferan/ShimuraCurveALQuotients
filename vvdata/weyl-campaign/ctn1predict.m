// N > 1: THE PREDICTIVE TEST.  Eight bases NEVER fitted, both supports each.
//
// The derivation (see theorem-Eeis-N.md) says: for D = D'*s squarefree with two
// primes and N a prime coprime to D,
//
//     rho_mu = 1/2 [Theta(D',Ns) - Theta(D',N)]/(mass_Ns - mass_N)
//              - 1/2 theta(DN,1)
//
// with the MULTIPLICATIVE mass, and the proof collapses it to
//     pred = c_{D'}^ram c_s^ram c_N^Eich   (N not| c),   0   (N | c)
// which is SYMMETRIC in D' <-> s.  So the prediction is twofold:
//   (i)  the identity holds at bases nobody has looked at, nothing fitted;
//   (ii) BOTH choices of D' give it -- support non-uniqueness at N > 1, the
//        exact analogue of what 1155/330/210 showed at N = 1.
//
// Control in the same run: the two-term part alone (third term dropped) must
// FAIL, so the test is not vacuous.
AttachSpec("ShimuraQuotients.spec");
load "vvdata/weyl-campaign/fastcosets.m";
load "vvdata/weyl-campaign/etared.m";
load "vvdata/weyl-campaign/ctlat.m";
load "vvdata/weyl-campaign/ctlatat.m";
PREC := 120; CC := ComplexField(PREC); QSIGN := -1;
frac := func< r | r - Floor(r) >;

massmul := function(Dp, R)
    ps := PrimeDivisors(Dp);
    m := Rationals()!(&*[ Integers() | p-1 : p in ps ]) / (48*2^(#ps-1));
    for p in PrimeDivisors(R) do m *:= Rationals()!(p+1)/2; end for;
    return m;
end function;

// D = p*q (two primes), N prime, N coprime to D.  NONE of these was ever fitted.
bases := [ <6,5>, <10,3>, <14,3>, <6,7>, <14,5>, <26,3>, <39,2>, <34,3> ];

nrun := 0; nhold := 0;
for b in bases do
    D := b[1]; N := b[2]; DN := D*N; lev := 4*DN;
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
    rho := ctThetaAt(Gm, words, CC, QSIGN, estc);
    rn := Maximum([ Abs(x) : x in rho ]);
    v3 := ctTheta(grossGram(DN, 1), words, CC, QSIGN);
    printf "\n%o_%o  DN=%o lev=%o  %o cosets  |A|=%o  #iso=%o  tracked order %o%o\n",
        D, N, DN, lev, nw, &*ms_, #iso, ord,
        ord eq N select " (= N, in the Eichler plane)" else " (NOT N -- check)";

    for Dp in PrimeDivisors(D) do
        s := D div Dp;
        Rl := N; Rls := N*s;
        m1 := massmul(Dp, Rl); msx := massmul(Dp, Rls);
        v1 := ctTheta(grossGram(Dp, Rl),  words, CC, QSIGN);
        vs := ctTheta(grossGram(Dp, Rls), words, CC, QSIGN);
        w1 := -(CC!m1)/(CC!(2*(msx-m1)));
        ws :=  (CC!msx)/(CC!(2*(msx-m1)));
        pred  := [ CC | ws*vs[i] + w1*v1[i] - (CC!(1/2))*v3[i] : i in [1..nw] ];
        ctrl  := [ CC | ws*vs[i] + w1*v1[i]                    : i in [1..nw] ];
        dev  := Maximum([ Abs(pred[i] - rho[i]) : i in [1..nw] ]);
        devc := Maximum([ Abs(ctrl[i] - rho[i]) : i in [1..nw] ]);
        ok := dev lt 10^(-90)*rn;
        nrun +:= 1; if ok then nhold +:= 1; end if;
        printf "   D'=%-3o s=%-3o  %o*T(%o,%o) %o*T(%o,%o) -1/2*T(%o,1)  worst %o  %o"
             cat "   [control, no third term: rel %o]\n",
            Dp, s, RationalField()!(-m1/(2*(msx-m1))), Dp, Rl,
            RationalField()!(msx/(2*(msx-m1))), Dp, Rls, DN,
            RealField(6)!dev, ok select "PREDICTED" else Sprintf("FAILS rel %o", RealField(6)!(dev/rn)),
            RealField(6)!(devc/rn);
    end for;
end for;
printf "\n%o of %o support-instances PREDICTED\n", nhold, nrun;
printf "DONE\n";
quit;
