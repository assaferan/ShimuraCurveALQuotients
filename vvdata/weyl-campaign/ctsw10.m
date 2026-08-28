// GAP 5d: p = 2 -- the 70 model bases the odd law cannot reach.
//
// Divide out everything now known (the global 12, the ramified factors at odd
// p | D, the Eichler factors at odd p | N) and read off what is left at 2:
//
//   R(m) = -a_E(m) / [ 12 H(4m) prod_{p|D, p odd} ram_p prod_{p|N, p odd} eich_p ]
//
// keyed by the 2-adic data of -4m = D_0 c^2, i.e. (v_2(c), (D_0|2)).  That is
// the parametrisation that worked at odd primes; m mod 8 was only ever a
// squarefree-m stand-in for it.
//
// Three regimes, because 2 plays three different roles and they must be
// separated:
//   2 | D, N = 1     ramified, no Eichler prime
//   2 | D, N odd     ramified, with an odd Eichler prime
//   D odd, N = 2     2 IS the Eichler prime
// The first two must agree if "ramified at 2" is one local factor.
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

sig := function(p, j)
    if j lt 0 then return Rationals()!0; end if;
    return Rationals()!((p^(j+1) - 1) div (p - 1));
end function;
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

//        D,  N,  s     regime label
cases := [ <10,1,2,"2|D, N=1">,  <14,1,2,"2|D, N=1">,  <22,1,2,"2|D, N=1">,
           <26,1,2,"2|D, N=1">,
           <10,3,2,"2|D, N odd">, <22,3,2,"2|D, N odd">, <10,7,2,"2|D, N odd">,
           <6,5,2,"2|D, N odd">,
           <15,2,3,"D odd, N=2">, <21,2,3,"D odd, N=2">, <35,2,5,"D odd, N=2">,
           <33,2,3,"D odd, N=2"> ];

GLOB := AssociativeArray();
for cs in cases do
    D := cs[1]; N := cs[2]; s := cs[3]; lab := cs[4]; Dp := D div s; DN := D*N;
    if N eq 1 then
        t1 := thetabar(Dp, 1); ts := thetabar(Dp, s);
        E := [ Rationals() | (((s+1)/2)*ts[i] - t1[i])/((s-1)/2) : i in [1..DEPTH] ];
    else
        m1 := massmul(Dp, N); msx := massmul(Dp, N*s);
        t1 := thetabar(Dp, N); ts := thetabar(Dp, N*s); t3 := thetabar(DN, 1);
        E := [ Rationals() | (msx/(2*(msx-m1)))*ts[i] - (m1/(2*(msx-m1)))*t1[i]
                             - (1/2)*t3[i] : i in [1..DEPTH] ];
    end if;

    R := AssociativeArray();
    for m in [1..DEPTH-1] do
        d0 := FundamentalDiscriminant(-4*m);
        c := Isqrt((-4*m) div d0);
        kn := Rationals()!12*Hur(4*m);
        for p in PrimeDivisors(D) do
            if p eq 2 then continue; end if;
            k := Valuation(c,p); chi := KroneckerSymbol(d0,p);
            kn *:= (1 - chi)/((p-1)*(sig(p,k) - chi*sig(p,k-1)));
        end for;
        for p in PrimeDivisors(N) do
            if p eq 2 then continue; end if;
            k := Valuation(c,p); chi := KroneckerSymbol(d0,p);
            kn *:= (p - chi)*p^k/((p^2-1)*(sig(p,k) - chi*sig(p,k-1)));
        end for;
        if kn eq 0 then
            if E[m+1] ne 0 then printf "  ** %o_%o m=%o: known 0 but a_E=%o\n", D,N,m,E[m+1]; end if;
            continue;
        end if;
        key := <Valuation(c,2), KroneckerSymbol(d0,2)>;
        if not IsDefined(R, key) then R[key] := {@ @}; end if;
        Include(~R[key], -E[m+1]/kn);
    end for;
    nc := #[ k : k in Keys(R) | #R[k] eq 1 ];
    printf "\n%o_%-3o [%o]  %o of %o (v2(c), (D0|2)) classes constant\n",
        D, N, lab, nc, #Keys(R);
    for key in Sort([ k : k in Keys(R) ]) do
        printf "   v2(c)=%o chi2=%2o : %o\n", key[1], key[2],
            #R[key] eq 1 select Sprintf("%o", Rep(R[key]))
            else Sprintf("%o distinct %o", #R[key],
                 #R[key] le 4 select [x : x in R[key]] else "(many)");
        if #R[key] eq 1 then
            gk := <lab, key[1], key[2]>;
            if not IsDefined(GLOB, gk) then GLOB[gk] := {@ @}; end if;
            Include(~GLOB[gk], Rep(R[key]));
        end if;
    end for;
end for;

printf "\n=== the 2-factor by (regime, v2(c), chi2), across bases ===\n";
for gk in Sort([ k : k in Keys(GLOB) ]) do
    printf "  %-12o v2=%o chi2=%2o : %o\n", gk[1], gk[2], gk[3],
        #GLOB[gk] eq 1 select Sprintf("%o", Rep(GLOB[gk]))
        else Sprintf("VARIES ACROSS BASES %o", [x : x in GLOB[gk]]);
end for;
printf "DONE\n";
quit;
