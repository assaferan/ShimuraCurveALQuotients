// GAP 5b: NON-SQUAREFREE m.  Read off the conductor corrections.
//
// Prop 9.15 covers m squarefree.  That is not enough to evaluate the pairing:
// sum_m c_f(-m) W(m) runs over EVERY m up to f's pole order, and those include
// 4, 8, 9, 12, ... .  So the formula is an insight, not yet an evaluator.
//
// Parametrise the way H does.  Write -4m = D_0 c^2 with D_0 FUNDAMENTAL.  For
// squarefree m, (D_0|p) = (-m|p) at odd p and c is 1 or 2, so the known law is
//
//   L(m) = prod_{p|D} (1 - (D_0|p))/(p-1) * 12(N - (D_0|N))/(N^2 - 1)
//
// written in terms of D_0.  Then divide out and look at what is left:
//
//   R(m) = -a_E(m) / ( H(4m) * L(m) )
//
// keyed by v_p(c) for p | DN.  R should be 1 where every v_p(c) = 0, and the
// deviations are the conductor corrections at the bad primes.  Bases are all-odd
// DN so that 2 is a GOOD prime and cannot confound the reading.
AttachSpec("ShimuraQuotients.spec");
load "vvdata/weyl-campaign/ctlat.m";
DEPTH := 150;

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

//        D,  N, D', R, Rs      -- DN odd throughout, so 2 is a good prime
cases := [ <15,7,5,7,21>, <35,3,7,3,15>, <21,5,7,5,15> ];

for cs in cases do
    D := cs[1]; N := cs[2]; Dp := cs[3]; Rl := cs[4]; Rls := cs[5]; DN := D*N;
    m1 := massmul(Dp, Rl); msx := massmul(Dp, Rls);
    t1 := thetabar(Dp, Rl); ts := thetabar(Dp, Rls); t3 := thetabar(DN, 1);
    w1 := -m1/(2*(msx-m1)); ws := msx/(2*(msx-m1));
    G := [ Rationals() | ws*ts[i] + w1*t1[i] - (1/2)*t3[i] : i in [1..DEPTH] ];
    ps := Sort(PrimeDivisors(DN));

    printf "\n=== %o_%o   D = %o, N = %o,  bad primes %o ===\n", D, N, D, N, ps;
    R := AssociativeArray();
    for m in [1..DEPTH-1] do
        d0 := FundamentalDiscriminant(-4*m);
        c := Isqrt((-4*m) div d0);
        L := (&*[ Rationals() | 1 - KroneckerSymbol(d0,p) : p in PrimeDivisors(D) ])
             / (&*[ Rationals() | p-1 : p in PrimeDivisors(D) ])
             * (Rationals()!12*(N - KroneckerSymbol(d0,N))/(N^2-1));
        if L eq 0 then
            if G[m+1] ne 0 then printf "  ** m=%o: L = 0 but a_E = %o\n", m, G[m+1]; end if;
            continue;
        end if;
        key := [ Valuation(c, p) : p in ps ];
        if not IsDefined(R, key) then R[key] := {@ @}; end if;
        Include(~R[key], -G[m+1]/(Hur(4*m)*L));
    end for;
    nc := 0;
    for key in Sort([ k : k in Keys(R) ]) do
        if #R[key] eq 1 then nc +:= 1; end if;
        printf "  v_p(c) = %-12o :  %o distinct  %o\n",
            Sprintf("%o", key), #R[key],
            #R[key] le 5 select [ x : x in R[key] ] else "(many)";
    end for;
    printf "  --> %o of %o conductor classes constant\n", nc, #Keys(R);
end for;
printf "DONE\n";
quit;
