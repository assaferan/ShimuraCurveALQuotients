// WHICH RATIONAL FIBRE OF THE 2030 STAR, AND WHICH TWIST (the 2002 stepH-enumerate, transposed).
// Y = X_0^2030(1)/<w_10, w_14, w_58> is hyperelliptic (certificate, PR #142); its branch divisor is
// a fibre of the star's hyperelliptic map; C has the rational Weierstrass point e = -3/2, so
// C : y^2 = q1(u) with q1 a quintic, u = (x - x0)/(x - e), and Y = Y_c : y^2 = q1(t^2/c),
// c supported on {-1, 2, 5, 7, 29}.  #Y_c(F_p) = #C(F_p) + #C'_c(F_p) - (p+1), C'_c : v^2 = c u q1(u).
AttachSpec("ShimuraQuotients.spec");
SetColumns(0);
BND := assigned BOUND select StringToInteger(BOUND) else 1000000;
S := {2, 5, 7, 29};
P<x> := PolynomialRing(Rationals());
C  := HyperellipticCurve(P![39,149,159,9,13,29,1], P![0,1,1]);
Cs := SimplifiedModel(C);  g := HyperellipticPolynomials(Cs);
e  := -3/2;  assert Evaluate(g, e) eq 0;
H := {1, 10, 14, 35, 58, 145, 203, 2030};
G := { m : m in Divisors(2030) | GCD(m, 2030 div m) eq 1 };
Yq := CreateShimuraQuot(2030, 1, H);  Yq`g := GenusShimuraCurveQuotient(2030, 1, H);
Xq := CreateShimuraQuot(2030, 1, G);  Xq`g := GenusShimuraCurveQuotient(2030, 1, G);
primes := [ p : p in PrimesUpTo(47) | 2030 mod p ne 0 ];
tY := [ ComputePointsViaTrace(Yq, p, 1) : p in primes ];
tX := [ ComputePointsViaTrace(Xq, p, 1) : p in primes ];
mX := [ #Points(ChangeRing(Cs, GF(p))) : p in primes ];
printf "star model vs trace at %o: %o\n", primes, tX eq mX;
need := [ tY[i] - tX[i] + primes[i] + 1 : i in [1..#primes] ];
printf "#Y(F_p), trace    %o\nrequired #C'(F_p) %o\n", tY, need;
pts := Points(Cs : Bound := BND);
fib := Sort(Setseq({ pt[1]/pt[3] : pt in pts | pt[3] ne 0 and pt[2] ne 0 }));
printf "%o rational points of height <= %o, %o finite non-Weierstrass fibres: %o\n", #pts, BND, #fib, fib;
PT<t> := PolynomialRing(Rationals());
twists := Sort([ s * &*[ Integers() | q : q in T ] : T in Subsets(S), s in [1, -1] ]);
hits := [];
for k in [0..#fib] do
    if k eq 0 then
        x0 := "infinity";
        yw := [ pt[2] : pt in pts | pt[3] eq 0 ][1];
        q1 := &+[ Coefficient(g, i) * (e*t + 1)^i * t^(6-i) : i in [0..6] ];
    else
        x0 := fib[k];
        yw := [ pt[2] : pt in pts | pt[3] ne 0 and pt[1]/pt[3] eq x0 ][1];
        q1 := &+[ Coefficient(g, i) * (e*t - x0)^i * (t - 1)^(6-i) : i in [0..6] ];
    end if;
    assert Degree(q1) eq 5 and Coefficient(q1, 0) ne 0;
    sunit := Set(PrimeDivisors(Integers() ! yw)) subset S;
    printf "\nfibre x = %-9o weighted y = %-8o S-unit: %o", x0, yw, sunit;
    for c in twists do
        Cp := ReducedMinimalWeierstrassModel(HyperellipticCurve(c * t * q1));
        ok := true;
        for i in [1..#primes] do
            try n := #Points(ChangeRing(Cp, GF(primes[i]))); catch err n := need[i]; end try;
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
