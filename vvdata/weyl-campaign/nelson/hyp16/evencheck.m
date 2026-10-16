// The sixteen all carry an Atkin-Lehner involution u whose quotient has genus 1 (genus 3) or 2
// (genus 4).  With the hyperelliptic involution iota central, u acts on the x-line as an
// involution, so x can be chosen with u : x -> -x and the octic is EVEN: f(x) = g(x^2).  Then
// X/u is the genus-1 curve Y^2 = g(t) and X/(iota u) the genus-2 curve w^2 = t g(t).
// Checked here against the stored models of X_0^35(1), where all three are known.
SetColumns(0);
P<t> := PolynomialRing(Rationals());
f0 := P![ -7/256, 0, -69/320, 0, -277/80000, 0, -1/40000, 0, -7/100000000 ];   // X_0^35(1)
g  := P![ -7/256, -69/320, -277/80000, -1/40000, -7/100000000 ];               // f0(x) = g(x^2)
e7 := P![ -7/256, -345/64, -277/128, -25/64, -7/256 ];                         // stored X/w_7
h5 := P![ 0, -7/256, -345/64, -277/128, -25/64, -7/256 ];                      // t*e7, genus 2?
printf "f0(x) = g(x^2) : %o\n", Evaluate(g, t^2) eq f0;
for lam in [1/25] do
  printf "g(t) equals the stored X_0^35(1)/w_7 quartic with t -> t/25 : %o\n",
     g eq Evaluate(e7, lam*t);
end for;
printf "stored X/w_7 genus %o ; w^2 = t*g(t) genus %o (stored X/w_5 genus 2)\n",
   Genus(HyperellipticCurve(e7)), Genus(HyperellipticCurve(t*g));
printf "X/w_5 stored:  %o\n", HyperellipticCurve(P![ 0, -7/256, -345/64, -277/128, -25/64, -7/256 ]);
C5 := HyperellipticCurve(P![ 0, -7/256, -345/64, -277/128, -25/64, -7/256 ]);
Cq := HyperellipticCurve(t*g);
printf "w^2 = t g(t) isomorphic over Q to the stored X_0^35(1)/w_5 (quintic form): %o\n",
   IsIsomorphic(Cq, C5);
// the three low-height invariants of the quartic with the two marked points t = 0, infinity
c := Coefficients(g);
A := c[1]*c[3]/c[2]^2; B := c[2]*c[4]/c[3]^2; C := c[3]*c[5]/c[4]^2;
printf "\nmarked-quartic invariants (unchanged by t -> lambda t and by scaling g):\n";
printf "   A = c0c2/c1^2 = %o\n   B = c1c3/c2^2 = %o\n   C = c2c4/c3^2 = %o\n", A, B, C;
printf "heights: %o %o %o\n", Maximum(Abs(Numerator(A)),Abs(Denominator(A))),
   Maximum(Abs(Numerator(B)),Abs(Denominator(B))), Maximum(Abs(Numerator(C)),Abs(Denominator(C)));
printf "EVENCHECK DONE\n";
quit;
