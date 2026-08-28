// N > 1: THE IDENTITY, with no fitted scalar.
//
// The residue lambda = 2N/(N+1) was NOT arithmetic.  The mass helper was called
// as mass([D'], [Rs]) with Rs the WHOLE Eichler level, so it contributed one
// factor (Rs+1)/2 instead of the product over the PRIMES of Rs.  At N = 1 the
// support is (D',1) and (D',s) with s PRIME, so the two agree and twenty bases
// never saw it.  At N > 1 the level Rs = R*s is composite and they diverge --
// and lambda was exactly the fudge that hid it.  It could only ever be a fudge:
// a single scalar cannot repair a class-dependent error, which is why it failed
// on {q,4,8,qr,4r,8r}.
//
// With the multiplicative mass the claim carries NO free parameter:
//
//   rho_mu = 1/2 * [Theta(D',Rs) - Theta(D',R)]/(mass_Rs - mass_R)
//            - 1/2 * theta(DN,1),      mass(D',R) = mass(D',1) * prod_{p|R} (p+1)/2
//
// and it reproduces the banked empirical weights exactly, e.g. at 22_3
//   1/2*(3 theta_6 - 2 theta_3) - 1/2 theta_66 = -theta(11,3) + 3/2 theta(11,6)
//                                                - 1/2 theta(66,1).
AttachSpec("ShimuraQuotients.spec");
load "vvdata/weyl-campaign/fastcosets.m";
load "vvdata/weyl-campaign/etared.m";
load "vvdata/weyl-campaign/ctlat.m";
load "vvdata/weyl-campaign/ctlatat.m";
PREC := 120; CC := ComplexField(PREC); QSIGN := -1;
frac := func< r | r - Floor(r) >;

// mass(D', R) with R factored into PRIMES -- the only change from the old helper
massmul := function(Dp, R)
    ps := PrimeDivisors(Dp);
    m := Rationals()!(&*[ Integers() | p-1 : p in ps ]) / (48*2^(#ps-1));
    for p in PrimeDivisors(R) do m *:= Rationals()!(p+1)/2; end for;
    return m;
end function;

//        D,  N, D', R, Rs
cases := [ <15,2,5,2,6>, <21,2,7,2,6>, <22,3,11,3,6>, <33,2,11,2,6>,
           <35,2,7,2,10>, <55,2,11,2,10>, <77,2,11,2,14>, <10,7,5,7,14> ];

nfail := 0;
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

    rho := ctThetaAt(Gm, words, CC, QSIGN, estc);
    m1 := massmul(Dp, Rl); msx := massmul(Dp, Rls);
    v1 := ctTheta(grossGram(Dp, Rl),  words, CC, QSIGN);
    vs := ctTheta(grossGram(Dp, Rls), words, CC, QSIGN);
    v3 := ctTheta(grossGram(DN, 1),   words, CC, QSIGN);
    w1 := -(CC!m1)/(CC!(2*(msx-m1)));            // weight on theta(D',R)
    ws :=  (CC!msx)/(CC!(2*(msx-m1)));           // weight on theta(D',Rs)
    w3 := CC!(-1/2);                             // weight on theta(DN,1)
    pred := [ CC | ws*vs[i] + w1*v1[i] + w3*v3[i] : i in [1..nw] ];

    rn  := Maximum([ Abs(x) : x in rho ]);
    dev := Maximum([ Abs(pred[i] - rho[i]) : i in [1..nw] ]);
    ok  := dev lt 10^(-90)*rn;
    if not ok then nfail +:= 1; end if;
    printf "%o_%-3o  weights  %o*T(%o,%o) %o*T(%o,%o) %o*T(%o,1)   %o cosets  worst %o  %o\n",
        D, N, RationalField()!(-m1/(2*(msx-m1))), Dp, Rl,
        RationalField()!(msx/(2*(msx-m1))), Dp, Rls,
        RationalField()!(-1/2), DN, nw, RealField(6)!dev,
        ok select "HOLDS" else Sprintf("MISMATCH rel %o", RealField(6)!(dev/rn));
    if not ok then
        badc := AssociativeArray();
        for i in [1..nw] do
            if Abs(pred[i] - rho[i]) lt 10^(-90)*rn then continue; end if;
            g := VVWordMatrix(words[i]); cl := GCD(g[2][1] mod lev, lev);
            badc[cl] := (IsDefined(badc,cl) select badc[cl] else 0) + 1;
        end for;
        printf "    disagree by class gcd(c,%o): %o\n", lev,
            [ <cl, badc[cl]> : cl in Sort([x : x in Keys(badc)]) ];
    end if;
end for;
printf "%o of %o bases HOLD\n", #cases - nfail, #cases;
printf "DONE\n";
quit;
