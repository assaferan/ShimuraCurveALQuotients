// N COMPOSITE probe: what do the isotropic cosets look like when the Shimura
// lattice is Eichler at TWO primes?  The N-prime proof needs mu to sit in the
// hyperbolic plane at N; with N = N1*N2 there is a 2-dimensional supply and the
// support of mu decides which combination the identity must match.
AttachSpec("ShimuraQuotients.spec");
frac := func< r | r - Floor(r) >;
for cs in [ <35,6>, <15,2> ] do
    D := cs[1]; N := cs[2];
    Ld := ShimuraCurveLattice(D, N);
    Gm := ChangeRing(Ld`Q, Rationals()); dn := Ld`denom;
    Gg := Ld`disc_grp; elts := [ g : g in Gg ];
    mods := Moduli(Gg); keep := [r : r in [1..#mods] | mods[r] gt 1];
    ms_ := [mods[r] : r in keep]; k := #ms_;
    Qb := [ frac(-(ChangeRing(e@@Ld`to_disc, Rationals())*Gm,
                   ChangeRing(e@@Ld`to_disc, Rationals()))/(2*dn^2)) : e in elts ];
    iso := [ i : i in [1..#elts] | Qb[i] eq 0 ];
    printf "\n=== D=%o N=%o  DN=%o  |A|=%o  moduli=%o  #iso=%o\n",
        D, N, D*N, &*ms_, ms_, #iso;
    for j in [1..#iso] do
        e := [ (Eltseq(elts[iso[j]])[keep[r]]) mod ms_[r] : r in [1..k] ];
        ord := LCM([ Integers() | ms_[r] div GCD(ms_[r], e[r]) : r in [1..k] ]);
        printf "   iso[%o] = %-16o order %-4o support %o%o\n", j, e, ord,
            PrimeDivisors(ord), j eq 2 select "   <- the convention's tracked coset" else "";
    end for;
end for;
printf "DONE\n";
quit;
