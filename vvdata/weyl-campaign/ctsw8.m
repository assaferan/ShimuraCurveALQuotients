// GAP 5b step 3: the formula for ALL m, tested predictively.
//
// Read off ctsw6.log with -4m = D_0 c^2, k_p = v_p(c), chi_p = (D_0|p), and
// sigma_j(p) = 1 + p + ... + p^j  (sigma_{-1} = 0):
//
//     conductor correction at p  =  p^{k [p = N]} / ( sigma_k - chi sigma_{k-1} )
//
// numerator 1 at a ramified prime, p^k at the Eichler prime; SAME denominator.
// At k = 0 both are 1, recovering Prop 9.15.  So for D odd, N an odd prime, and
// ANY m:
//
//   -a_E(m) = H(4m)
//             * prod_{p|D} (1 - chi_p)/((p-1) (sigma_{k_p} - chi_p sigma_{k_p-1}))
//             * 12 (N - chi_N) N^{k_N} / ((N^2-1)(sigma_{k_N} - chi_N sigma_{k_N-1}))
//
// Fitted on 15_7, 35_3, 21_5 with p in {3,5,7} and k <= 2.  Tested below on
// fresh N (11, 13), a fresh D, and to a depth that reaches k = 3.
AttachSpec("ShimuraQuotients.spec");
load "vvdata/weyl-campaign/ctlat.m";
DEPTH := 800;

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

sig := function(p, j)   // 1 + p + ... + p^j, and 0 for j < 0
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

//        D,  N, D', R, Rs        fitted?
cases := [ <15,7,5,7,21, true> ];

for cs in cases do
    D := cs[1]; N := cs[2]; Dp := cs[3]; Rl := cs[4]; Rls := cs[5];
    fitted := cs[6]; DN := D*N;
    m1 := massmul(Dp, Rl); msx := massmul(Dp, Rls);
    t1 := thetabar(Dp, Rl); ts := thetabar(Dp, Rls); t3 := thetabar(DN, 1);
    w1 := -m1/(2*(msx-m1)); ws := msx/(2*(msx-m1));
    G := [ Rationals() | ws*ts[i] + w1*t1[i] - (1/2)*t3[i] : i in [1..DEPTH] ];

    nok := 0; nbad := 0; bl := [ ]; kmax := 0; nsf := 0;
    for m in [1..DEPTH-1] do
        d0 := FundamentalDiscriminant(-4*m);
        c := Isqrt((-4*m) div d0);
        pr := Rationals()!1;
        for p in PrimeDivisors(D) do
            k := Valuation(c, p); chi := KroneckerSymbol(d0, p);
            kmax := Maximum(kmax, k);
            pr *:= (1 - chi)/((p-1)*(sig(p,k) - chi*sig(p,k-1)));
        end for;
        kN := Valuation(c, N); chiN := KroneckerSymbol(d0, N);
        kmax := Maximum(kmax, kN);
        pr *:= Rationals()!12*(N - chiN)*N^kN
               / ((N^2-1)*(sig(N,kN) - chiN*sig(N,kN-1)));
        pred := -Hur(4*m)*pr;
        if not IsSquarefree(m) then nsf +:= 1; end if;
        if G[m+1] eq pred then nok +:= 1;
        else nbad +:= 1; if #bl lt 4 then Append(~bl, <m, G[m+1], pred>); end if;
        end if;
    end for;
    printf "%o_%-3o %-9o %o of %o m < %o agree (%o of them non-squarefree, k up to %o)  %o\n",
        D, N, fitted select "(fitted)" else "(FRESH)", nok, nok+nbad, DEPTH, nsf, kmax,
        nbad eq 0 select "PREDICTED" else Sprintf("%o FAIL e.g. %o", nbad, bl);
end for;
printf "DONE\n";
quit;
