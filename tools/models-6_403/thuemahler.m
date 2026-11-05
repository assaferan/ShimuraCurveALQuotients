// THE BRANCH FIBRE OF Y = X_0^6(403)/<w_2,w_13,w_93> -> X_0^6(403)*, EXACTLY (the 2002 argument).
// Y has good reduction at p not dividing 6*403 = 2418; its branch points on the star C are distinct modulo
// such p; the branch divisor is a fibre {(a : +-y' : b)} of y^2 = g(x) (hyperellipticity is proved by
// the certificate of PR #142), and the two points agree modulo p exactly when p | y'.  So y' is a
// unit outside S = {2, 3, 13, 31}, hence so is G(a, b) = b^6 g(a/b) = y'^2 = (2a + 3b) Q5(a, b), and
// Res(x - 1, Q5(x, 1)) = 2^6 5 7^2 makes Q5(a, b) an S-unit: a Thue-Mahler equation of degree five,
// solved by Gherga-Siksek's solver (17 solutions, listed verbatim in thuemahler.sols, 246 s).
AttachSpec("ShimuraQuotients.spec");
SetColumns(0);
S := {2, 3, 13, 31};
sols := eval Read("tools/models-6_403/thuemahler.sols");
printf "%o listed solutions of the Thue-Mahler equation\n", #sols;
P<x> := PolynomialRing(Rationals());
C  := HyperellipticCurve(P![9,-21,33,-47,15,9,1], P![0,1,1]);
Cs := SimplifiedModel(C);  g := HyperellipticPolynomials(Cs);
printf "y^2 = g(x), g = %o;  primes of its discriminant %o\n", g, PrimeDivisors(Integers() ! Discriminant(Cs));
G  := func< a, b | &+[ Integers() ! Coefficient(g, i) * a^i * b^(6-i) : i in [0..6] ] >;
Q5 := func< a, b | 4*a^5 + 40*a^4*b + 101*a^3*b^2 - 85*a^2*b^3 + 48*a*b^4 - 36*b^5 >;
L  := func< a, b | a - b >;
assert forall{ <a, b> : a in [-3..3], b in [-3..3] | G(a, b) eq L(a, b) * Q5(a, b) };
printf "Res(x - 1, Q5(x, 1)) = %o\n", Factorisation(Integers() ! Resultant(P![-1,1], P![-36,48,-85,101,40,4]));
fibres := [];
for s in sols do
    a := s[1]; b := s[2];
    assert GCD(a, b) eq 1;
    assert Abs(Q5(a, b)) eq 2^s[3] * 3^s[4] * 13^s[5] * 31^s[6];
    sq, yw := IsSquare(G(a, b));
    if sq and yw ne 0 and Set(PrimeDivisors(Integers() ! yw)) subset S then Append(~fibres, <a, b, yw>); end if;
end for;
printf "%o of them are points of C with S-unit weighted y-coordinate (the possible branch fibres):\n   %o\n",
       #fibres, [ <b eq 0 select "oo" else Sprint(a/b), yw> where a := f[1] where b := f[2] where yw := f[3] : f in fibres ];
e := 1;
H := {1, 2, 13, 26, 93, 186, 1209, 2418};
Yq := CreateShimuraQuot(6, 403, H);  Yq`g := GenusShimuraCurveQuotient(6, 403, H);
primes := [ p : p in PrimesUpTo(47) | 2418 mod p ne 0 ];
need := [ ComputePointsViaTrace(Yq, p, 1) - #Points(ChangeRing(Cs, GF(p))) + p + 1 : p in primes ];
PT<t> := PolynomialRing(Rationals());
twists := Sort([ sg * &*[ Integers() | q : q in T ] : T in Subsets(S), sg in [1, -1] ]);
hits := [];
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
        if ok then Append(~hits, <a, b, c>); printf "MATCH: fibre (%o : %o), twist c = %o\n", a, b, c; end if;
    end for;
end for;
printf "%o matching (fibre, twist) pair(s) over %o fibres and %o twists\n", #hits, #fibres, #twists;
assert #hits eq 1 and hits[1][1] eq -7*hits[1][2] and hits[1][3] eq 1;   // the fibre x = -7 with twist 1
print "so the branch fibre is x = -7 and Y is y^2 = q1(t^2), q1 the quintic of the star at u = (x + 7)/(x - 1)";
quit;
