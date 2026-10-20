// Reconstructing a genus-3 quotient from an elliptic quotient E and the branch divisor alone
// (ExactDoubleCover.m), on quotients whose models are already stored in data/models:
//
//   X_0^35(1) and X_0^51(1): hyperelliptic, with the involution x -> -x of the stored even octic
//   y^2 = F(x) (the equations of Guo-Yang, Compositio 153 (2017), as stored in
//   data/models); E is y^2 = F_even(u) with u = x^2 and the branch divisor is the fibre over u = 0
//   and u = oo. The expected outcome is a single cover over Q, hyperelliptic, whose model after
//   the twist read off the trace-formula counts is Q-isomorphic to the stored octic. The same
//   cover written as a double cover of its canonical conic (the presentation that exists even when
//   the conic has no rational point) has the point counts the function field of the cover gives,
//   and so does a twist of that cover whose conic has no rational point, where no model y^2 = f(x)
//   over Q exists.
//
//   X_0^6(23)/w_23: not hyperelliptic. Its stored model was built in this repository as the fibre
//   product of its two elliptic quotients y1^2 = q1(x), y2^2 = q2(x) over the star, so this case
//   checks agreement with the trace formula and with that construction, not with a published
//   equation. Over E: y1^2 = q1 the branch divisor is the fibre over the two roots of q2 that are
//   not roots of q1. Here E(Q)[2] has order 2 and there are TWO covers over Q with this branch
//   divisor, one hyperelliptic and one not; the trace formula is consistent with the
//   non-hyperelliptic one only, and for it kappa*f/q2(x) is a square in the function field of E,
//   so the selected cover is the stored curve exactly; for the other cover no twist works.
//
// The trace-formula counts come from ComputePointsViaTrace.

P<x> := PolynomialRing(Rationals());
integral := func<H | HyperellipticCurve(P!(HyperellipticPolynomials(H) *
    LCM([Denominator(c) : c in Coefficients(HyperellipticPolynomials(H))])^2))>;

// An elliptic curve isomorphic to the genus-one curve y^2 = q(x), the isomorphism, and the
// function x carried over to it.
function elliptic_with_x(q)
    den := LCM([Denominator(c) : c in Coefficients(q)]);
    C := HyperellipticCurve(q*den^2);
    pts := Points(C : Bound := 5000);
    assert #pts gt 0;
    E, mp := EllipticCurve(C, pts[1]);
    return E, Pullback(Inverse(mp), FunctionField(C).1);
end function;

function trace_hits(E, covers, D, N, W)
    X := CreateShimuraQuot(D, N, W);
    X`g := GenusShimuraCurveQuotient(D, N, W);
    primes := [p : p in PrimesInInterval(3, 130) | (D*N) mod p ne 0];
    target := [ComputePointsViaTrace(X, p, 1) : p in primes];
    return MatchCoverToCounts(E, covers, target, primes : ExtraPrimes := PrimeDivisors(D*N));
end function;

for D in [35, 51] do
    models := eval (Read("data/models/models_" cat IntegerToString(D) cat "_1.m") cat " return models;");
    Ftop := models[[Integers()|1]][1][2];
    assert &and[Coefficient(Ftop, i) eq 0 : i in [1,3,5,7]];
    E, u := elliptic_with_x(P![Coefficient(Ftop, 2*i) : i in [0..4]]);
    Z, Pl := SignDecomposition(Divisor(u));
    covers := DoubleCoversWithBranchDivisor(E, Z + Pl);
    assert #covers eq 1 and covers[1]`genus eq 3 and covers[1]`noether eq 5;
    hits := trace_hits(E, covers, D, 1, {1});
    assert #hits eq 1;
    H := HyperellipticModelOfCover(TwistedCover(covers[1], hits[1][2]));
    assert IsIsomorphic(H, integral(HyperellipticCurve(Ftop)));
    // the same cover written over its canonical conic: where that model has smooth reduction its
    // point count must be the count of the curve, which the function field of the cover gives
    // independently, and PointCountOverConic must refuse at the other primes rather than answer
    sel := TwistedCover(covers[1], hits[1][2]);
    conic, G := HyperellipticModelOverConic(sel);
    pr := [p : p in PrimesInInterval(3, 60) | D mod p ne 0];
    viaff := CoverPointCounts(E, sel, pr);
    viaconic := [PointCountOverConic(conic, G, p) : p in pr];
    nusable := 0;
    for i in [1..#pr] do
        if viaff[i][1] lt 0 or viaconic[i] lt 0 then continue; end if;
        assert viaconic[i] eq viaff[i][1];
        nusable +:= 1;
    end for;
    assert nusable ge 8;
    printf "X_0^%o(1): one cover, hyperelliptic, twist %o, Q-isomorphic to the stored model\n", D, hits[1][2];
end for;

models := eval (Read("data/models/models_6_23.m") cat " return models;");
crv := models[[Integers()|1,23]][1];
R<s,z,y1,y2> := PolynomialRing(Rationals(), 4);
eqs := [R | eval e : e in crv[3]];
q1 := -Evaluate(eqs[1], [x, 1, 0, 0]);
q2 := -Evaluate(eqs[2], [x, 1, 0, 0]);
assert Degree(GCD(q1, q2)) eq 2;
E, xE := elliptic_with_x(q1);
Z := SignDecomposition(Divisor(Evaluate(q2 div GCD(q1, q2), xE)));
covers := DoubleCoversWithBranchDivisor(E, Z);
assert #covers eq 2 and {c`noether : c in covers} eq {5, 6} and &and[c`genus eq 3 : c in covers];
hits := trace_hits(E, covers, 6, 23, {1,23});
assert #hits eq 1 and covers[hits[1][1]]`noether eq 6;
// the selected cover, twisted by kappa, is the stored curve y2^2 = q2(x) over E: kappa*f/q2(x) is a
// square in the function field of E, and for the other cover no twist makes it one
sel := covers[hits[1][1]]; other := [c : c in covers | c`noether eq 5][1];
assert IsSquare(hits[1][2] * sel`f / Evaluate(q2, xE));
assert not exists{k : k in [1, -1, 2, -2, 3, -3, 6, -6, 23, -23, 46, -46, 69, -69, 138, -138] | IsSquare(k * other`f / Evaluate(q2, xE))};
printf "X_0^6(23)/w_23: two covers, the trace formula keeps the non-hyperelliptic one (twist %o), which is the stored curve\n", hits[1][2];

// A cover whose canonical conic has NO rational point: there is no model y^2 = f(x) over Q, and the
// model over the conic is the only presentation. The configuration is two points of E(Q(sqrt -2)) on
// 1122i1, the elliptic quotient of X_0^6(187)/<w_3, w_11, w_17>; it was found by the search over
// Heegner heights, but nothing here depends on that -- the test only needs a branch divisor whose
// cover has this shape. The count over the conic is compared with the count from the function field
// of the cover, which is computed by a different route and does not use the conic at all.
printf "A cover whose canonical conic has no rational point...";
E := EllipticCurve("1122i1");
K<r> := QuadraticField(-2);
B := Divisor(Place(E(K)![65/4 + 3*r, -131/8 + 165/8*r])) + Divisor(Place(E(K)![18 - 4*r, -6 - 28*r]));
covers := DoubleCoversWithBranchDivisor(E, B);
sel := [TwistedCover(c, -14) : c in covers | c`noether eq 5];
assert #sel eq 1;
conic, G, ws, yy := HyperellipticModelOverConic(sel[1]);
P2 := ProjectiveSpace(Rationals(), 2);
assert not HasRationalPoint(Conic(P2, Evaluate(conic, [P2.1, P2.2, P2.3])));
refused := false;
try _ := HyperellipticModelOfCover(sel[1]); catch err refused := true; end try;
assert refused;
pr := [p : p in PrimesInInterval(3, 60) | p ne 11];
cff := CoverPointCounts(E, sel[1], pr);
nchk := 0;
for i in [1..#pr] do
    cc := PointCountOverConic(conic, G, pr[i]);
    if cff[i][1] lt 0 or cc lt 0 then continue; end if;
    assert cc eq cff[i][1];
    nchk +:= 1;
end for;
assert nchk ge 6;
printf " the model over the conic exists and agrees with the function field at %o primes\n", nchk;
