// WHICH RATIONAL FIBRE OF THE STAR, AND WHICH TWIST  (paper, Section 3, parts (3) and (4)).
//
//   magma -b tools/proof-2002/stepH-enumerate.m < /dev/null          (about three minutes)
//   magma -b BOUND:=100000 tools/proof-2002/stepH-enumerate.m < /dev/null
//
// If Y = X_0^2002(1)/<w_2, w_13, w_77> is hyperelliptic, its branch divisor on the star C is a
// fibre {u = 0} of C's hyperelliptic map, C has a model y^2 = q1(u) with q1 a quintic, and Y is
// one of the twists Y_c : y^2 = q1(t^2/c), c supported on {-1, 2, 7, 11, 13}.  The Klein four-group
// <wbar, iota> on Y_c has quotients C, P^1 and C'_c : v^2 = c u q1(u), so
//         #Y_c(F_p) = #C(F_p) + #C'_c(F_p) - (p+1),
// and the twists are told apart by the genus-two curves C'_c.  The trace formula gives #Y(F_p).
//
// The script runs this over every fibre of C above a rational point of height at most BOUND, not
// only the fibre the numerics point to, and reports every (fibre, c) whose counts agree with Y's at
// the eleven primes.  It also prints the S-unit test of Remark 3 for each fibre: the two branch
// points cannot coincide modulo a prime of good reduction of Y, so the weighted y-coordinate of the
// branch fibre is a unit outside {2, 7, 11, 13}.

AttachSpec("ShimuraQuotients.spec");
SetColumns(0);
BND := assigned BOUND select StringToInteger(BOUND) else 1000000;
S := {2, 7, 11, 13};

P<x> := PolynomialRing(Rationals());
C  := HyperellipticCurve(P![ 4, 64, 196, -456, 180, 4, 16 ], P![ 0, 1, 1 ]);   // the proved star
Cs := SimplifiedModel(C);                                                       // y^2 = 4f + h^2
g  := HyperellipticPolynomials(Cs);
e  := -1/8;  assert Evaluate(g, e) eq 0;                                        // the rational Weierstrass point
printf "4f + h^2 = %o\n   factors as %o\n", g, Factorisation(g);
printf "Jac(C)(Q)[2] has order %o\n", #TwoTorsionSubgroup(Jacobian(Cs));

H := {1, 2, 13, 26, 77, 154, 1001, 2002};
G := { m : m in Divisors(2002) | GCD(m, 2002 div m) eq 1 };
Yq := CreateShimuraQuot(2002, 1, H);  Yq`g := GenusShimuraCurveQuotient(2002, 1, H);
Xq := CreateShimuraQuot(2002, 1, G);  Xq`g := GenusShimuraCurveQuotient(2002, 1, G);
primes := [ p : p in PrimesUpTo(47) | 2002 mod p ne 0 ];
tY := [ ComputePointsViaTrace(Yq, p, 1) : p in primes ];
tX := [ ComputePointsViaTrace(Xq, p, 1) : p in primes ];
mX := [ #Points(ChangeRing(Cs, GF(p))) : p in primes ];
assert tX eq mX;     // the star model and the trace formula agree, as Theorem 1 says they must
need := [ tY[i] - tX[i] + primes[i] + 1 : i in [1..#primes] ];
printf "\nprimes            %o\n#Y(F_p), trace    %o\nrequired #C'(F_p) %o\n", primes, tY, need;

pts := Points(Cs : Bound := BND);
fib := Sort(Setseq({ pt[1]/pt[3] : pt in pts | pt[3] ne 0 and pt[2] ne 0 }));
printf "\n%o rational points of C of height <= %o, in %o finite non-Weierstrass fibres and the fibre at infinity\n",
       #pts, BND, #fib;

PT<t> := PolynomialRing(Rationals());
twists := Sort([ s * &*[ Integers() | q : q in T ] : T in Subsets(S), s in [1, -1] ]);
hits := [];
for k in [0..#fib] do
    if k eq 0 then
        x0 := "infinity";
        yw := [ pt[2] : pt in pts | pt[3] eq 0 ][1];
        q1 := &+[ Coefficient(g, i) * (e*t + 1)^i * t^(6-i) : i in [0..6] ];     // u = 1/(x - e)
    else
        x0 := fib[k];
        yw := [ pt[2] : pt in pts | pt[3] ne 0 and pt[1]/pt[3] eq x0 ][1];
        q1 := &+[ Coefficient(g, i) * (e*t - x0)^i * (t - 1)^(6-i) : i in [0..6] ]; // u = (x - x0)/(x - e)
    end if;
    assert Degree(q1) eq 5 and Coefficient(q1, 0) ne 0;
    sunit := Set(PrimeDivisors(Integers() ! yw)) subset S;
    printf "\nfibre x = %-9o weighted y = %-8o S-unit: %o", x0, yw, sunit;
    for c in twists do
        Cp := ReducedMinimalWeierstrassModel(HyperellipticCurve(c * t * q1));
        ok := true;
        for i in [1..#primes] do
            try n := #Points(ChangeRing(Cp, GF(primes[i]))); catch err n := need[i]; end try;  // singular: no information
            if n ne need[i] then ok := false; break; end if;
        end for;
        if ok then
            Append(~hits, <x0, c>);
            printf "\n   MATCH with c = %o:  C' = %o\n   Y_c : y^2 = %o", c, Cp, Evaluate(q1, t^2/c);
        end if;
    end for;
end for;
printf "\n\n%o matching (fibre, twist) pair(s) over %o fibres and %o twists\n", #hits, #fib + 1, #twists;
quit;
