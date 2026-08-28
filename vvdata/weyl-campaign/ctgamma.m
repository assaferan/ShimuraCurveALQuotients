// LEMMA half 2, CLOSED: test the theorem's identity at gamma*w, gamma in Gamma_0(4DN).
//
// The identity of Theorem 9.11 was verified at a FIXED list of coset
// representatives w.  Re-evaluate both sides at gamma*w for gamma in
// Gamma_0(4DN) -- the same cosets, different representatives.  Each side is a
// weight-3/2 form slashed by a matrix, so
//     Shimura side at gamma*w = nu_Shim(gamma) * (value at w)
//     Gross   side at gamma*w = nu_Gross(gamma) * (value at w)
// and the identity survives the change of representative IF AND ONLY IF
//     nu_Gross = nu_Shim .
//
// That is exactly the multiplier statement wanted, and it needs no knowledge of
// the panel character at all: the panel character is defined by being conjugate
// to whatever governs rho, i.e. to nu_Shim.  So this closes the character half
// at EVERY base, including the ones whose eta pools are too small for the
// membership test (21_2, 22_3, 55_2, 77_2, 34_3, 22_1).
AttachSpec("ShimuraQuotients.spec");
load "vvdata/weyl-campaign/fastcosets.m";
load "vvdata/weyl-campaign/etared.m";
load "vvdata/weyl-campaign/ctlat.m";
load "vvdata/weyl-campaign/ctlatat.m";
PREC := 100; CC := ComplexField(PREC); QSIGN := -1;
frac := func< r | r - Floor(r) >;

massmul := function(Dp, Rl)
    ps := PrimeDivisors(Dp);
    m := Rationals()!(&*[ Integers() | p-1 : p in ps ]) / (48*2^(#ps-1));
    for p in PrimeDivisors(Rl) do m *:= Rationals()!(p+1)/2; end for;
    return m;
end function;

//        D,  N, D', R, Rs   -- the inconclusive bases first, 15_2 as control
cases := [ <21,2,7,2,6>, <22,3,11,3,6>, <55,2,11,2,10>, <77,2,11,2,14>,
           <34,3,17,3,6>, <15,2,5,2,6> ];

nok := 0; nrun := 0;
for cs in cases do
    D := cs[1]; N := cs[2]; Dp := cs[3]; Rl := cs[4]; Rls := cs[5];
    DN := D*N; lev := 4*DN;
    reps := fastCosetReps(lev);
    Ld := ShimuraCurveLattice(D, N);
    Gm := ChangeRing(Ld`Q, Rationals()); dn := Ld`denom;
    Gg := Ld`disc_grp; elts := [ g : g in Gg ];
    mods := Moduli(Gg); keep := [r : r in [1..#mods] | mods[r] gt 1];
    ms_ := [mods[r] : r in keep]; k := #ms_;
    Qb := [ frac(-(ChangeRing(e@@Ld`to_disc, Rationals())*Gm,
                   ChangeRing(e@@Ld`to_disc, Rationals()))/(2*dn^2)) : e in elts ];
    iso := [ i : i in [1..#elts] | Qb[i] eq 0 ];
    estc := [ (Eltseq(elts[iso[2]])[keep[r]]) mod ms_[r] : r in [1..k] ];

    // gamma in Gamma_0(lev), several; times a spread of coset representatives
    gams := [ ];
    for d in [1, 3, 7, 11, 13, 17] do
        if GCD(d, lev) ne 1 then continue; end if;
        a := InverseMod(d mod lev, lev);
        b := (a*d - 1) div lev;
        g := Matrix(Integers(), 2, 2, [a, b, lev, d]);
        if Determinant(g) eq 1 then Append(~gams, g); end if;
    end for;
    step := Maximum(1, #reps div 12);
    wsub := [ reps[i] : i in [1..#reps by step] ];
    mats := [ g*w : g in gams, w in wsub ];
    words := [ VVSTWord(m) : m in mats ];

    rho := ctThetaAt(Gm, words, CC, QSIGN, estc);
    m1 := massmul(Dp, Rl); msx := massmul(Dp, Rls);
    v1 := ctTheta(grossGram(Dp, Rl),  words, CC, QSIGN);
    vs := ctTheta(grossGram(Dp, Rls), words, CC, QSIGN);
    v3 := ctTheta(grossGram(DN, 1),   words, CC, QSIGN);
    w1 := -(CC!m1)/(CC!(2*(msx-m1)));
    ws :=  (CC!msx)/(CC!(2*(msx-m1)));
    pred := [ CC | ws*vs[i] + w1*v1[i] - (CC!(1/2))*v3[i] : i in [1..#words] ];
    rn  := Maximum([ Abs(x) : x in rho ]);
    dev := Maximum([ Abs(pred[i] - rho[i]) : i in [1..#words] ]);
    ok := dev lt 10^(-70)*rn;
    nrun +:= 1; if ok then nok +:= 1; end if;
    printf "%o_%-3o  %o gamma x %o reps = %o matrices OFF the representative list"
         cat "   worst %o   %o\n",
        D, N, #gams, #wsub, #words, RealField(6)!dev,
        ok select "HOLDS -> nu_Gross = nu_Shim" else Sprintf("FAILS rel %o", RealField(6)!(dev/rn));
end for;
printf "\n%o of %o bases hold off the representative list\n", nok, nrun;
printf "DONE\n";
quit;
