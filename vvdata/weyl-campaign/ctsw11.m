// GAP 5 COMPLETE: ONE formula, every prime including 2, every m, real bases.
//
// ctsw10 read off the 2-adic factors and they are NOT a special case -- they are
// the odd law with p = 2 substituted in:
//     ram_p  = (1 - chi_p) / ((p-1)(sigma_k - chi sigma_{k-1}))
//     eich_p = (p - chi_p) p^k / ((p^2-1)(sigma_k - chi sigma_{k-1}))
// e.g. ram at p=2: k=1 chi=-1 -> 2/4 = 1/2, k=2 chi=0 -> 1/7, k=3 chi=-1 ->
// 2/22 = 1/11; eich at p=2: k=0 chi=0 -> 2/3, k=1 chi=0 -> 4/9, k=2 chi=1 ->
// 4/12 = 1/3.  All matched.
//
// The earlier "p = 2 is an exception" was an artefact of keying on m mod 8 while
// restricted to squarefree m.  Under -4m = D_0 c^2 the tables
// delta_2 = 1, 1/2, 0 and 8, 8, 6, 8, 8, 4 are just this law's shadow.
//
// So the claim is now: for D squarefree with omega(D) = 2, N = 1 or prime
// coprime to D (2 allowed in EITHER slot), and ANY m >= 1,
//
//   -a_E(m) = 12 H(4m) prod_{p|D} ram_p prod_{p|N} eich_p .
//
// Tested below on real model bases across all three regimes.
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
printf "Hurwitz gate: passed\n\n";

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

//     D, N, s   -- all real entries of data/models except 15_7 (marked)
cases := [ <15,1,3>, <35,1,5>, <39,1,3>, <55,1,5>,          // D odd,  N=1
           <10,1,2>, <14,1,2>, <22,1,2>, <26,1,2>,          // 2|D,    N=1
           <15,2,3>, <21,2,3>, <33,2,3>, <35,2,5>,          // D odd,  N=2
           <10,3,2>, <22,3,2>, <10,7,2>, <6,5,2>,           // 2|D,    N odd
           <14,5,2>, <34,3,2>,
           <15,7,3> ];                                      // synthetic control

tot := 0; totbad := 0;
for cs in cases do
    D := cs[1]; N := cs[2]; s := cs[3]; Dp := D div s; DN := D*N;
    if N eq 1 then
        t1 := thetabar(Dp, 1); ts := thetabar(Dp, s);
        E := [ Rationals() | (((s+1)/2)*ts[i] - t1[i])/((s-1)/2) : i in [1..DEPTH] ];
    else
        m1 := massmul(Dp, N); msx := massmul(Dp, N*s);
        t1 := thetabar(Dp, N); ts := thetabar(Dp, N*s); t3 := thetabar(DN, 1);
        E := [ Rationals() | (msx/(2*(msx-m1)))*ts[i] - (m1/(2*(msx-m1)))*t1[i]
                             - (1/2)*t3[i] : i in [1..DEPTH] ];
    end if;
    nok := 0; nbad := 0; bl := [ ];
    for m in [1..DEPTH-1] do
        d0 := FundamentalDiscriminant(-4*m);
        c := Isqrt((-4*m) div d0);
        pr := Rationals()!12*Hur(4*m);
        for p in PrimeDivisors(D) do
            k := Valuation(c,p); chi := KroneckerSymbol(d0,p);
            pr *:= (1 - chi)/((p-1)*(sig(p,k) - chi*sig(p,k-1)));
        end for;
        for p in PrimeDivisors(N) do
            k := Valuation(c,p); chi := KroneckerSymbol(d0,p);
            pr *:= (p - chi)*p^k/((p^2-1)*(sig(p,k) - chi*sig(p,k-1)));
        end for;
        if E[m+1] eq -pr then nok +:= 1;
        else nbad +:= 1; if #bl lt 3 then Append(~bl, <m, E[m+1], -pr>); end if;
        end if;
    end for;
    tot +:= nok + nbad; totbad +:= nbad;
    printf "%3o_%-3o (s=%-2o) %-14o %o of %o  %o\n", D, N, s,
        IsEven(D) select (N eq 1 select "[2|D, N=1]" else "[2|D, N odd]")
                  else (N eq 1 select "[D odd, N=1]" else (N eq 2 select "[D odd, N=2]" else "[all odd]")),
        nok, nok+nbad,
        nbad eq 0 select "ALL m AGREE" else Sprintf("%o FAIL e.g. %o", nbad, bl);
end for;
printf "\n%o of %o (base, m) pairs agree across %o bases\n", tot-totbad, tot, #cases;
printf "DONE\n";
quit;
