// omega(D) = 4 at N > 1: the last scope gap of theorem-Eeis-N.md.
// D = 210 (four ramified primes), N = 11, level 9240 -- the smallest instance.
// Two choices of s, the first two nonzero isotropic cosets each.
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

//         D,  N     (N composite; 22_3 is the N-prime regression check)
bases := [ <210,11> ];

ninst := 0; nhold := 0;
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
    printf "\n===== D=%o N=%o  DN=%o  lev=%o  %o cosets  #iso=%o =====\n",
        D, N, DN, lev, nw, #iso;

    mus := [ ]; Ss := [ ];
    for j in [2..#iso] do
        e := [ (Eltseq(elts[iso[j]])[keep[r]]) mod ms_[r] : r in [1..k] ];
        ord := LCM([ Integers() | ms_[r] div GCD(ms_[r], e[r]) : r in [1..k] ]);
        Append(~mus, e); Append(~Ss, PrimeDivisors(ord));
    end for;

    for s in PrimeDivisors(D)[1..2] do
        Dp := D div s;
        printf "  --- s = %o, D' = %o ---\n", s, Dp;
        // Term(T) for EVERY T subset of the primes of N, precomputed.  (Magma
        // will not let a function assign to a variable of the enclosing scope,
        // so there is no closure-cached version of this.)
        tcache := AssociativeArray();
        for T in Subsets(SequenceToSet(PrimeDivisors(N))) do
            Ts := Sort(SetToSequence(T));
            key := &*[ Integers() | p : p in Ts ];
            R := N div key;
            if IsOdd(#Ts) then
                v := ctTheta(grossGram(Dp*s*key, R), words, CC, QSIGN);
            else
                m1 := massmul(Dp*key, R); m2 := massmul(Dp*key, R*s);
                v1 := ctTheta(grossGram(Dp*key, R),   words, CC, QSIGN);
                v2 := ctTheta(grossGram(Dp*key, R*s), words, CC, QSIGN);
                v := [ CC | (CC!m2*v2[i] - CC!m1*v1[i])/CC!(m2-m1) : i in [1..nw] ];
            end if;
            tcache[key] := v;
        end for;

        for j in [1..Minimum(2,#mus)] do
            S := Ss[j];
            rho := ctThetaAt(Gm, words, CC, QSIGN, mus[j]);
            rn := Maximum([ Abs(x) : x in rho ]);
            pred := [ CC | 0 : i in [1..nw] ];
            desc := [ ];
            for T in Subsets(SequenceToSet(S)) do
                Ts := Sort(SetToSequence(T));
                v := tcache[&*[ Integers() | p : p in Ts ]];
                sg := IsOdd(#Ts) select -1 else 1;
                w := CC!sg / CC!(2^#S);
                pred := [ CC | pred[i] + w*v[i] : i in [1..nw] ];
                key := &*[ Integers() | p : p in Ts ];
                Append(~desc, IsOdd(#Ts)
                    select Sprintf("%o/%o T(%o,%o)", sg, 2^#S, Dp*s*key, N div key)
                    else Sprintf("%o/%o Tel(%o;%o,%o)", sg, 2^#S, Dp*key, N div key, (N div key)*s));
            end for;
            dev := Maximum([ Abs(pred[i] - rho[i]) : i in [1..nw] ]);
            ok := dev lt 10^(-90)*rn;
            ninst +:= 1; if ok then nhold +:= 1; end if;
            printf "   mu=%-16o S=%-10o %o   worst %o  %o\n",
                mus[j], S, desc, RealField(6)!dev,
                ok select "HOLDS" else Sprintf("MISMATCH rel %o", RealField(6)!(dev/rn));
        end for;
    end for;
end for;
printf "\n%o of %o (base, s, mu) instances HOLD\n", nhold, ninst;
printf "DONE\n";
quit;
