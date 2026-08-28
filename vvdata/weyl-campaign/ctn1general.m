// THE GENERAL N > 1 IDENTITY -- N squarefree with ANY number of primes, and
// EVERY nonzero isotropic tracked coset, not just the convention's iso[2].
//
// Let D be squarefree with omega(D) even, D = D'*s with s a prime; let N be
// squarefree, coprime to D; let mu be a nonzero isotropic coset of A_L and
// S = { p | N : mu_p != 0 }.  Then
//
//     ct_L(gamma, mu) = 2^(-|S|) * sum_{T subset S} (-1)^|T| Term(T)
//
//     Term(T) = theta( D'*s*prod T , N/prod T )                    |T| ODD
//             = Tel over s of (D'*prod T, N/prod T) and
//                             (D'*prod T, (N/prod T)*s)            |T| EVEN
//
// (|T| even makes omega(D'*s*prod T) even, i.e. INDEFINITE and unavailable as a
// Gross lattice, which is exactly when the telescoping over s is needed.)
//
// WHY.  The Weil entry factors prime by prime with a lattice-independent global
// factor (ctsplit.m).  mu isotropic gives kappa = prod_{p in S} [p not| c] and
// NO phase.  And
//     prod_{p in S} (c_p^Eich - c_p^ram)/2 = prod_{p in S} c_p^Eich * kappa ,
// since the difference is 2/p off the bad cosets and 0 on them.  Expanding that
// product over subsets T is precisely the alternating sum above.
//
// The proven N-prime case is S = {N}, T in {{}, {N}}: Tel((D',N),(D',Ns)) minus
// theta(D's N, 1) = theta(DN,1), over 2 -- the banked three-term identity.
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
bases := [ <22,3>, <35,6>, <15,14> ];

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

    for s in PrimeDivisors(D) do
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

        for j in [1..#mus] do
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
