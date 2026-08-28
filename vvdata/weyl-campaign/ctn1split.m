// N > 1: read off the LOCAL factors of all four lattices, prime by prime.
// Self-check first: glob * prod_p loc_p must reproduce ctTheta exactly.
AttachSpec("ShimuraQuotients.spec");
load "vvdata/weyl-campaign/fastcosets.m";
load "vvdata/weyl-campaign/etared.m";
load "vvdata/weyl-campaign/ctlat.m";
load "vvdata/weyl-campaign/ctsplit.m";
PREC := 120; CC := ComplexField(PREC); QSIGN := -1;
RF := RealField(PREC); tol := RF!10^(-80);
frac := func< r | r - Floor(r) >;
R8 := ComplexField(8);

//        D,  N, D', R, Rs
cases := [ <22,3,11,3,6>, <15,2,5,2,6>, <10,7,5,7,14> ];

for cs in cases do
    D := cs[1]; N := cs[2]; Dp := cs[3]; Rl := cs[4]; Rls := cs[5];
    DN := D*N; lev := 4*DN;
    words := [ VVSTWord(g) : g in fastCosetReps(lev) ];
    nw := #words;
    Ld := ShimuraCurveLattice(D, N);
    Gm := ChangeRing(Ld`Q, Rationals());
    printf "\n=== %o_%o   DN=%o  lev=%o ===\n", D, N, DN, lev;

    lats := [* <"L(shimura)", Gm>, <Sprintf("G(%o,%o)", Dp, Rl), grossGram(Dp, Rl)>,
               <Sprintf("G(%o,%o)", Dp, Rls), grossGram(Dp, Rls)>,
               <Sprintf("G(%o,1)", DN), grossGram(DN, 1)> *];
    nl := #lats;
    SP := [* *]; CT := [* *];
    for j in [1..nl] do
        okj, globj, locj, plistj, msj, sigj := ctSplitLocal(lats[j][2], words, CC, QSIGN);
        ctj := ctTheta(lats[j][2], words, CC, QSIGN);
        // self-check
        bad := 0; worst := RF!0;
        for i in [1..nw] do
            if not okj[i] then continue; end if;
            pr := globj[i];
            for p in plistj do pr *:= locj[i][p]; end for;
            d := Abs(pr - ctj[i]);
            if d gt tol then bad +:= 1; end if;
            if d gt worst then worst := d; end if;
        end for;
        printf "  %-12o moduli=%o  sig8=%o  self-check: %o bad, worst %o\n",
            lats[j][1], msj, [ <p, sigj[p]> : p in plistj ],
            bad, RealField(6)!worst;
        Append(~SP, <okj, globj, locj, plistj>);
        Append(~CT, ctj);
    end for;

    // one representative per cusp class, among the words the split covers
    reps := AssociativeArray();
    for i in [1..nw] do
        if not SP[1][1][i] then continue; end if;
        g := VVWordMatrix(words[i]); c := g[2][1];
        gg := GCD(c mod lev, lev);
        if not IsDefined(reps, gg) then reps[gg] := i; end if;
    end for;
    pall := Sort(SetToSequence(SequenceToSet(&cat[ SP[j][4] : j in [1..nl] ])));
    printf "  --- loc_p by class (one representative each) ---\n";
    printf "  %-6o %-14o %o\n", "class", "lattice", [ Sprintf("p=%o", p) : p in pall ];
    for gg in Sort([x : x in Keys(reps)]) do
        i := reps[gg];
        for j in [1..nl] do
            row := [ ];
            for p in pall do
                Append(~row, p in SP[j][4] select Sprintf("%o", R8!SP[j][3][i][p])
                                             else "  --  ");
            end for;
            printf "  %-6o %-14o %o\n", j eq 1 select gg else "", lats[j][1], row;
        end for;
    end for;
end for;
printf "DONE\n";
quit;
