// THE BRANCH FIBRE, EXACTLY: every rational fibre of the star that the cover can be branched along,
// from a Thue-Mahler equation  (paper, Section 3, Proposition on the branch fibre).
//
//   magma -b tools/proof-2002/stepH-thuemahler.m < /dev/null         (under a minute)
//
// Y = X_0^2002(1)/<w_2, w_13, w_77> has good reduction at every p not dividing 2002, and its
// quotient map to the star C reduces modulo such p to a tame double cover of smooth curves, whose two
// branch points are therefore distinct modulo p.  If Y is hyperelliptic the branch divisor is a
// fibre {(a : +-y' : b)} of C's hyperelliptic map, in weighted coordinates on y^2 = g(x) with
// gcd(a, b) = 1, and the two points agree modulo p exactly when p | y'.  So y' is a unit outside
// S = {2, 7, 11, 13}, hence so is G(a, b) = b^6 g(a/b) = y'^2; and G = (8a + b) Q5(a, b) with the
// resultant of the two factors equal to 2^7 7 11 13, so Q5(a, b) is an S-unit: a Thue-Mahler
// equation of degree five.  Its solutions were computed by Gherga-Siksek's solver and are listed in
// stepH-thuemahler.sols; this script checks each listed pair is a solution, keeps those that are
// points of C with S-unit y', and runs the twist-and-count test of stepH-enumerate on each.  What it
// cannot check is the completeness of the list, which is the solver's theorem.

AttachSpec("ShimuraQuotients.spec");
SetColumns(0);
S := {2, 7, 11, 13};
sols := eval Read("tools/proof-2002/stepH-thuemahler.sols");
printf "%o listed solutions of the Thue-Mahler equation\n", #sols;

P<x> := PolynomialRing(Rationals());
C  := HyperellipticCurve(P![ 4, 64, 196, -456, 180, 4, 16 ], P![ 0, 1, 1 ]);
Cs := SimplifiedModel(C);  g := HyperellipticPolynomials(Cs);
printf "y^2 = g(x), g = %o;  primes of its discriminant %o\n", g, PrimeDivisors(Integers() ! Discriminant(Cs));
G  := func< a, b | &+[ Integers() ! Coefficient(g, i) * a^i * b^(6-i) : i in [0..6] ] >;
Q5 := func< a, b | 8*a^5 + a^4*b + 90*a^3*b^2 - 239*a^2*b^3 + 128*a*b^4 + 16*b^5 >;
L  := func< a, b | 8*a + b >;
assert forall{ <a, b> : a in [-3..3], b in [-3..3] | G(a, b) eq L(a, b) * Q5(a, b) };
Lx := P ! [1, 8]; Qx := P ! [16, 128, -239, 90, 1, 8];
printf "Res(8x + 1, Q5(x, 1)) = %o\n", Factorisation(Integers() ! Resultant(Lx, Qx));

fibres := [];
for s in sols do
    a := s[1]; b := s[2];
    assert GCD(a, b) eq 1;
    assert Abs(Q5(a, b)) eq 2^s[3] * 7^s[4] * 11^s[5] * 13^s[6];      // it is a solution
    sq, yw := IsSquare(G(a, b));
    if sq and yw ne 0 and Set(PrimeDivisors(Integers() ! yw)) subset S then
        Append(~fibres, <a, b, yw>);
    end if;
end for;
printf "%o of them are points of C with S-unit weighted y-coordinate (the possible branch fibres):\n   %o\n",
       #fibres, [ <b eq 0 select "oo" else Sprint(a/b), yw> where a := f[1] where b := f[2] where yw := f[3] : f in fibres ];

// the twist-and-count test along each of them
e := -1/8;
H := {1, 2, 13, 26, 77, 154, 1001, 2002};
Yq := CreateShimuraQuot(2002, 1, H);  Yq`g := GenusShimuraCurveQuotient(2002, 1, H);
primes := [ p : p in PrimesUpTo(47) | 2002 mod p ne 0 ];
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
assert hits eq [ <1, 0, 2> ];
print "so if Y is hyperelliptic, its branch fibre is the fibre at infinity and Y is Y_2";
quit;
