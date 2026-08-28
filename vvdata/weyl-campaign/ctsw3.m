// GAP 5, step 4: divide out the part that is now known, and READ OFF the rest.
//
// Step 3 (ctsw2.log) refuted the 15_2-fitted form at every other base -- but the
// failures were UNIFORM ratios per base: 2/3 at 21_2, 2/5 at 33_2, 1/3 at 35_2,
// which are 8/(2*6), 8/(2*10), 8/(4*6).  So the local factor at an odd ramified
// prime is (1 - (-m|p))/(p-1), not (1 - (-m|p)).  (The cases that "agreed" at
// those bases were the ones where both sides vanish.)
//
// Rather than guess the remaining factors at 2 and at the Eichler prime, divide
// out what is known,
//
//     R(m) = a_E(m) / [ -H(4m) * prod_{p | D, p odd} (1 - (-m|p))/(p-1) ] ,
//
// and print R keyed by the local data at 2 and at N.  Bases are chosen to
// separate the variables: D odd with N = 2, D odd with N odd, and 2 | D.
AttachSpec("ShimuraQuotients.spec");
load "vvdata/weyl-campaign/ctlat.m";
DEPTH := 160;

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

//        D,  N, D', R, Rs                             regime
cases := [ <15,2,5,2,6>,   <21,2,7,2,6>,               // D odd, N = 2
           <33,2,11,2,6>,  <35,2,7,2,10>,
           <15,7,5,7,21>,  <35,3,7,3,15>,              // D odd, N odd
           <22,3,11,3,6>,  <10,7,5,7,14> ];            // 2 | D

for cs in cases do
    D := cs[1]; N := cs[2]; Dp := cs[3]; Rl := cs[4]; Rls := cs[5]; DN := D*N;
    m1 := massmul(Dp, Rl); msx := massmul(Dp, Rls);
    t1 := thetabar(Dp, Rl); ts := thetabar(Dp, Rls); t3 := thetabar(DN, 1);
    w1 := -m1/(2*(msx-m1)); ws := msx/(2*(msx-m1));
    G := [ Rationals() | ws*ts[i] + w1*t1[i] - (1/2)*t3[i] : i in [1..DEPTH] ];
    oddD := [ p : p in PrimeDivisors(D) | IsOdd(p) ];

    printf "\n=== %o_%o   D = %o = %o,  Eichler N = %o%o\n", D, N, D,
        PrimeDivisors(D), N, IsEven(D) select "   (2 | D)" else "";
    R := AssociativeArray();
    for m in [1..DEPTH-1] do
        if not IsSquarefree(m) then continue; end if;
        den := Hur(4*m) * (&*[ Rationals() | 1 - KroneckerSymbol(-m,p) : p in oddD ]
                           / (&*[ Rationals() | p-1 : p in oddD ]));
        if den eq 0 then
            if G[m+1] ne 0 then printf "   ** m = %o: known part 0 but a_E = %o\n", m, G[m+1]; end if;
            continue;
        end if;
        key := <m mod 8, KroneckerSymbol(-m, N)>;
        if not IsDefined(R, key) then R[key] := {@ @}; end if;
        Include(~R[key], -G[m+1]/den);
    end for;
    nc := 0;
    for key in Sort([ k : k in Keys(R) ]) do
        if #R[key] eq 1 then nc +:= 1; end if;
        printf "   m=%o mod 8, (-m|%o)=%2o :  %o  %o\n", key[1], N, key[2], #R[key],
            #R[key] le 4 select [ x : x in R[key] ] else "(many)";
    end for;
    printf "   --> %o of %o classes constant\n", nc, #Keys(R);
end for;
printf "DONE\n";
quit;
