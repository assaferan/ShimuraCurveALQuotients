// GAP 5b step 2: key on BOTH v_p(c) and (D_0|p), and read off the conductor
// correction as a function of (p, k, chi).
//
// ctsw5 keyed only on v_p(c) and got R = 1 on every class with c coprime to DN
// -- the squarefree law recovered exactly -- but 2 or 3 values where p | c.
// That is the missing coordinate: once p | c the symbol (D_0|p) is still free,
// and the local factor depends on both.
AttachSpec("ShimuraQuotients.spec");
load "vvdata/weyl-campaign/ctlat.m";
DEPTH := 250;

Hur := function(n)
    if n eq 0 then return Rationals()!(-1)/12; end if;
    if n mod 4 in {1, 2} then return Rationals()!0; end if;
    s := Rationals()!0; f := 1;
    while f^2 le n do
        if n mod f^2 eq 0 then
            d := -(n div f^2);
            if d mod 4 in {0, 1} then
                w := d eq -3 select 6 else (d eq -4 select 4 else 2);
                s +:= Rationals()!(2*ClassNumber(QuadraticForms(d)))/w;
            end if;
        end if;
        f +:= 1;
    end while;
    return s;
end function;
T3 := ThetaSeries(StandardLattice(3), 4*DEPTH);
error if exists{ n : n in [1..90] |
    (n mod 4 in {1,2} and Coefficient(T3,n) ne 12*Hur(4*n)) or
    (n mod 8 eq 3 and Coefficient(T3,n) ne 24*Hur(n)) }, "Hurwitz routine wrong";
printf "Hurwitz gate: passed\n";

massmul := function(Dp, Rl)
    ps := PrimeDivisors(Dp);
    m := Rationals()!(&*[ Integers() | p-1 : p in ps ]) / (48*2^(#ps-1));
    for p in PrimeDivisors(Rl) do m *:= Rationals()!(p+1)/2; end for;
    return m;
end function;
thetabar := function(Dp, Rl)
    A := ChangeRing(grossGram(Dp, Rl), Rationals());
    L := LatticeWithGram(A);
    reps := Representatives(Genus(L));
    raw := [ Rationals() | 0 : m in [0..2*DEPTH] ];
    den := Rationals()!0;
    for Lr in reps do
        w := Rationals()!1/#AutomorphismGroup(Lr);
        T := ThetaSeries(Lr, 2*DEPTH);
        for m in [0..2*DEPTH] do raw[m+1] +:= w*Coefficient(T, m); end for;
        den +:= w;
    end for;
    raw := [ x/den : x in raw ];
    return [ raw[2*m+1] : m in [0..DEPTH-1] ];
end function;

cases := [ <15,7,5,7,21>, <35,3,7,3,15>, <21,5,7,5,15> ];

// accumulate per (role, p, k, chi) across bases, to see if the factor is local
GLOB := AssociativeArray();
for cs in cases do
    D := cs[1]; N := cs[2]; Dp := cs[3]; Rl := cs[4]; Rls := cs[5]; DN := D*N;
    m1 := massmul(Dp, Rl); msx := massmul(Dp, Rls);
    t1 := thetabar(Dp, Rl); ts := thetabar(Dp, Rls); t3 := thetabar(DN, 1);
    w1 := -m1/(2*(msx-m1)); ws := msx/(2*(msx-m1));
    G := [ Rationals() | ws*ts[i] + w1*t1[i] - (1/2)*t3[i] : i in [1..DEPTH] ];
    ps := Sort(PrimeDivisors(DN));

    printf "\n=== %o_%o   ramified %o,  Eichler %o ===\n", D, N, PrimeDivisors(D), N;
    R := AssociativeArray();
    for m in [1..DEPTH-1] do
        d0 := FundamentalDiscriminant(-4*m);
        c := Isqrt((-4*m) div d0);
        L := (&*[ Rationals() | 1 - KroneckerSymbol(d0,p) : p in PrimeDivisors(D) ])
             / (&*[ Rationals() | p-1 : p in PrimeDivisors(D) ])
             * (Rationals()!12*(N - KroneckerSymbol(d0,N))/(N^2-1));
        if L eq 0 then continue; end if;
        key := [ <p, Valuation(c,p), KroneckerSymbol(d0,p)> : p in ps ];
        if not IsDefined(R, key) then R[key] := {@ @}; end if;
        Include(~R[key], -G[m+1]/(Hur(4*m)*L));
    end for;
    nc := 0;
    for key in Sort([ k : k in Keys(R) ]) do
        if #R[key] eq 1 then nc +:= 1; end if;
        // only print classes where some prime actually divides c
        if exists{ t : t in key | t[2] gt 0 } then
            printf "  %-34o : %o  %o\n",
                Sprintf("%o", [ <t[1],t[2],t[3]> : t in key | t[2] gt 0 ]),
                #R[key], #R[key] le 4 select [ x : x in R[key] ] else "(many)";
        end if;
        // record single-bad-prime classes for the cross-base comparison
        act := [ t : t in key | t[2] gt 0 ];
        if #act eq 1 and #R[key] eq 1 then
            p := act[1][1];
            role := p eq N select "Eich" else "ram";
            gk := <role, p, act[1][2], act[1][3]>;
            if not IsDefined(GLOB, gk) then GLOB[gk] := {@ @}; end if;
            Include(~GLOB[gk], Rep(R[key]));
        end if;
    end for;
    printf "  --> %o of %o classes constant\n", nc, #Keys(R);
end for;

printf "\n=== the correction factor by (role, p, k = v_p(c), chi = (D_0|p)) ===\n";
for gk in Sort([ k : k in Keys(GLOB) ]) do
    printf "  %-5o p=%-3o k=%o chi=%2o :  %o\n", gk[1], gk[2], gk[3], gk[4],
        #GLOB[gk] eq 1 select Sprintf("%o", Rep(GLOB[gk]))
        else Sprintf("INCONSISTENT ACROSS BASES %o", [x : x in GLOB[gk]]);
end for;
printf "DONE\n";
quit;
