// THE BRANCH FIBRE OF Y = X_0^6(341)/<w_3,w_11,w_62> -> X_0^6(341)*, EXACTLY, GIVEN THE STAR MODEL.
// In the coordinate x -> 1/x the star is y^2 = g(x) with g = x (x - 3) Q4(x), Q4 = 36x^4 - 24x^3 - 35x^2
// + 13x - 8 irreducible, S = {2, 3, 11, 31}; the resultants of the factors are S-units, so at a fibre
// with S-unit weighted y-coordinate Q4(a, b) is an S-unit (both signs: degree 4 is even).  The two
// rational Weierstrass points x = 0, 3 give two rational theta characteristics, so for each fibre
// BOTH normal forms u = (x - x0)/(x - e), e in {0, 3}, must be tried.
AttachSpec("ShimuraQuotients.spec");
SetColumns(0);
S := {2, 3, 11, 31};
sols := eval Read("tools/models-6_341/thuemahler.sols");
printf "%o listed solutions of the Thue-Mahler equation (both signs)\n", #sols;
P<x> := PolynomialRing(Rationals());
C  := HyperellipticCurve(P![9,-33,9,29,-12,6], P![0,1,1]);
g  := P![0,24,-47,118,37,-132,36];
assert IsIsomorphic(HyperellipticCurve(g), C);
Cs := HyperellipticCurve(g);
G  := func< a, b | &+[ Integers() ! Coefficient(g, i) * a^i * b^(6-i) : i in [0..6] ] >;
Q4 := func< a, b | 36*a^4 - 24*a^3*b - 35*a^2*b^2 + 13*a*b^3 - 8*b^4 >;
assert forall{ <a, b> : a in [-3..3], b in [-3..3] | G(a, b) eq a*(a - 3*b)*Q4(a, b) };
fibres := [];
for s in sols do
    a := s[1]; b := s[2];
    assert GCD(a, b) eq 1;
    assert Abs(Q4(a, b)) eq 2^s[3] * 3^s[4] * 11^s[5] * 31^s[6];
    sq, yw := IsSquare(G(a, b));
    // (a, b) and (-a, -b) are the same fibre (Q4 is even in (a, b)); keep one representative
    if sq and yw ne 0 and Set(PrimeDivisors(Integers() ! yw)) subset S
       and not exists{ f : f in fibres | f[1]*b eq f[2]*a } then Append(~fibres, <a, b, yw>); end if;
end for;
printf "%o of them are points of C with S-unit weighted y-coordinate:\n   %o\n",
       #fibres, [ <b eq 0 select "oo" else Sprint(a/b), yw> where a := f[1] where b := f[2] where yw := f[3] : f in fibres ];
H := {1, 3, 11, 33, 62, 186, 682, 2046};
Yq := CreateShimuraQuot(6, 341, H);  Yq`g := GenusShimuraCurveQuotient(6, 341, H);
primes := [ p : p in PrimesUpTo(47) | 2046 mod p ne 0 ];
need := [ ComputePointsViaTrace(Yq, p, 1) - #Points(ChangeRing(Cs, GF(p))) + p + 1 : p in primes ];
PT<t> := PolynomialRing(Rationals());
twists := Sort([ sg * &*[ Integers() | q : q in T ] : T in Subsets(S), sg in [1, -1] ]);
hits := [];
for e in [0, 3] do
for f in fibres do
    a := f[1]; b := f[2];
    if b eq 0 then
        q1 := &+[ Coefficient(g, i) * (e*t + 1)^i * t^(6-i) : i in [0..6] ];
    else
        q1 := &+[ Coefficient(g, i) * (e*t - a/b)^i * (t - 1)^(6-i) : i in [0..6] ];
    end if;
    assert Degree(q1) eq 5 and Coefficient(q1, 0) ne 0;
    for c in twists do
        Cp := ReducedMinimalWeierstrassModel(HyperellipticCurve(c * t * q1));
        ok := true;
        for i in [1..#primes] do
            try n := #Points(ChangeRing(Cp, GF(primes[i]))); catch err n := need[i]; end try;
            if n ne need[i] then ok := false; break; end if;
        end for;
        if ok then Append(~hits, <e, a, b, c>); printf "MATCH: base e = %o, fibre (%o : %o), twist c = %o\n", e, a, b, c; end if;
    end for;
end for;
end for;
printf "%o matching (base, fibre, twist) triple(s) over 2 bases, %o fibres and %o twists\n", #hits, #fibres, #twists;
assert #hits eq 1 and hits[1][1] eq 3 and hits[1][2]*3 eq hits[1][3] and hits[1][4] eq -3;   // the fibre x = 1/3 with base e = 3 and twist -3
Y := HyperellipticCurve(P![81, 0, 468, 0, 1126, 0, 1420, 0, -47, 0, -744]);
q1 := &+[ Coefficient(g, i) * (3*x - 1/3)^i * (x - 1)^(6-i) : i in [0..6] ];
printf "y^2 = q1(t^2/-3) is y^2 = -744x^10 - 47x^8 + 1420x^6 + 1126x^4 + 468x^2 + 81 over Q: %o\n",
       IsIsomorphic(HyperellipticCurve(Evaluate(q1, x^2/-3)), Y);
print "so, given the star model, the branch fibre is x = 1/3 (x = 3 on the star as written) and Y is the curve above";
quit;
