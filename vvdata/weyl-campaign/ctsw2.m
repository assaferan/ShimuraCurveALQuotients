// GAP 5, step 3: the closed form, tested as a PREDICTION at fresh bases.
//
// Fitted at 15_2 only (ctsw15b.log, 24 of 24 local classes constant):
//
//     a_E(m) = - H(4m) * lambda_2(m) * prod_{p | D} ( 1 - (-m|p) )     m squarefree
//
// with lambda_2(m) = 1 (m = 1,2,5,6 mod 8), 3/4 (m = 3 mod 8), 1/2 (m = 7 mod 8).
// The odd factor 1 - (-m|p) is 2 at inert p, 1 at p | m, and 0 at SPLIT p, so
// the support rule -- a_E(m) vanishes unless every ramified prime of D is
// non-split in Q(sqrt(-m)), i.e. unless that field embeds in B_D -- is part of
// the formula rather than an extra condition.
//
// Nothing is refitted below.  The bases here have different D and, at 10_7 and
// 22_3, a different Eichler prime and 2 | D instead of 2 | N -- the case the
// 15_2 fit cannot see.
AttachSpec("ShimuraQuotients.spec");
load "vvdata/weyl-campaign/ctlat.m";
DEPTH := 200;

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

lam2 := function(m)
    case m mod 8:
        when 3: return Rationals()!3/4;
        when 7: return Rationals()!1/2;
        else return Rationals()!1;
    end case;
end function;

//        D,  N, D', R, Rs
cases := [ <15,2,5,2,6>, <21,2,7,2,6>, <33,2,11,2,6>, <35,2,7,2,10>,
           <22,3,11,3,6>, <10,7,5,7,14> ];

for cs in cases do
    D := cs[1]; N := cs[2]; Dp := cs[3]; Rl := cs[4]; Rls := cs[5]; DN := D*N;
    m1 := massmul(Dp, Rl); msx := massmul(Dp, Rls);
    t1 := thetabar(Dp, Rl); ts := thetabar(Dp, Rls); t3 := thetabar(DN, 1);
    w1 := -m1/(2*(msx-m1)); ws := msx/(2*(msx-m1));
    G := [ Rationals() | ws*ts[i] + w1*t1[i] - (1/2)*t3[i] : i in [1..DEPTH] ];

    nok := 0; nbad := 0; badlist := [ ];
    for m in [1..DEPTH-1] do
        if not IsSquarefree(m) then continue; end if;
        pred := -Hur(4*m)*lam2(m)*(&*[ Rationals() |
                    1 - KroneckerSymbol(-m, p) : p in PrimeDivisors(D) ]);
        if G[m+1] eq pred then nok +:= 1;
        else nbad +:= 1; if #badlist lt 6 then Append(~badlist, <m, G[m+1], pred>); end if;
        end if;
    end for;
    printf "%o_%-3o D = %-4o (primes %o), Eichler %o :  %o of %o squarefree m agree%o\n",
        D, N, D, PrimeDivisors(D), N, nok, nok + nbad,
        nbad eq 0 select "   PREDICTED"
        else Sprintf("   %o FAIL, e.g. <m, actual, pred> = %o", nbad, badlist);
end for;
printf "DONE\n";
quit;
