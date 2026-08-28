// GAP D: what the tracked-coset choice does at COMPOSITE N.
//
// The record says mult = 1/2 c_eta(0) for ANY nonzero isotropic eta, "they all
// agree", equivalently sum_iso c_eta(0) / (4(N-1)) with 2N-2 nonzero isotropic
// cosets.  Both counts are N-PRIME counts.  The general theorem says
//
//     rho^S = ct_L(gamma, 0) * prod_{p in S} [ p not| c ] ,   S = supp(mu),
//
// because prod_{p in S}(c_p^Eich - c_p^ram)/2 = prod_{p in S} c_p^Eich [p not|c]
// and the remaining factors are exactly ct_L(gamma,0).  So the tracked cosets do
// not give unrelated vectors: they give ONE vector restricted to different sets
// of cosets, nested by S.  This checks that, counts the classes, and reports how
// many DISTINCT functionals a base actually has.
AttachSpec("ShimuraQuotients.spec");
load "vvdata/weyl-campaign/fastcosets.m";
load "vvdata/weyl-campaign/etared.m";
load "vvdata/weyl-campaign/ctlat.m";
load "vvdata/weyl-campaign/ctlatat.m";
PREC := 80; CC := ComplexField(PREC); QSIGN := -1;
RF := RealField(PREC); tol := RF!10^(-60);
frac := func< r | r - Floor(r) >;

bases := [ <22,3>, <35,6>, <15,14> ];   // prime-N control, then two composite

for b in bases do
    D := b[1]; N := b[2]; DN := D*N; lev := 4*DN;
    words := [ VVSTWord(g) : g in fastCosetReps(lev) ];
    nw := #words;
    cvals := [ VVWordMatrix(w)[2][1] : w in words ];
    Ld := ShimuraCurveLattice(D, N);
    Gm := ChangeRing(Ld`Q, Rationals()); dn := Ld`denom;
    Gg := Ld`disc_grp; elts := [ g : g in Gg ];
    mods := Moduli(Gg); keep := [r : r in [1..#mods] | mods[r] gt 1];
    ms_ := [mods[r] : r in keep]; k := #ms_;
    Qb := [ frac(-(ChangeRing(e@@Ld`to_disc, Rationals())*Gm,
                   ChangeRing(e@@Ld`to_disc, Rationals()))/(2*dn^2)) : e in elts ];
    iso := [ i : i in [1..#elts] | Qb[i] eq 0 ];
    rho0 := ctTheta(Gm, words, CC, QSIGN);

    printf "\n=== %o_%o  N = %o = %o   #iso = %o (prod (2p-1) = %o)\n",
        D, N, N, PrimeDivisors(N), #iso,
        &*[ Integers() | 2*p-1 : p in PrimeDivisors(N) ];

    bySupp := AssociativeArray(); nbad := 0;
    funcs := [* *];   // distinct rho-vectors actually seen
    for j in [2..#iso] do
        e := [ (Eltseq(elts[iso[j]])[keep[r]]) mod ms_[r] : r in [1..k] ];
        ord := LCM([ Integers() | ms_[r] div GCD(ms_[r], e[r]) : r in [1..k] ]);
        S := PrimeDivisors(ord);
        key := &*[ Integers() | p : p in S ];
        bySupp[key] := (IsDefined(bySupp, key) select bySupp[key] else 0) + 1;
        rmu := ctThetaAt(Gm, words, CC, QSIGN, e);
        // predicted: rho0 restricted to the cosets with p not| c for all p in S
        ok := true;
        for i in [1..nw] do
            keepit := forall{ p : p in S | cvals[i] mod p ne 0 };
            pred := keepit select rho0[i] else CC!0;
            if Abs(rmu[i] - pred) gt tol then ok := false; end if;
        end for;
        if not ok then nbad +:= 1; end if;
        // is this vector new?
        isnew := true;
        for v in funcs do
            if forall{ i : i in [1..nw] | Abs(v[i] - rmu[i]) le tol } then
                isnew := false; break;
            end if;
        end for;
        if isnew then Append(~funcs, rmu); end if;
    end for;
    printf "  rho^S = rho^0 restricted to {p not| c, p in S}: %o of %o cosets agree\n",
        #iso - 1 - nbad, #iso - 1;
    printf "  by support: %o   (predicted prod_{p in S}(2p-2): %o)\n",
        [ <key, bySupp[key]> : key in Sort([x : x in Keys(bySupp)]) ],
        [ <key, &*[ Integers() | 2*p-2 : p in PrimeDivisors(key) ]>
          : key in Sort([x : x in Keys(bySupp)]) ];
    printf "  DISTINCT rho-vectors among the %o nonzero isotropic cosets: %o"
         cat "   (2^omega(N) - 1 = %o)\n",
        #iso - 1, #funcs, 2^(#PrimeDivisors(N)) - 1;
    printf "  --> \"they all agree\" is %o here\n",
        #funcs eq 1 select "TRUE" else "FALSE";
end for;
printf "DONE\n";
quit;
