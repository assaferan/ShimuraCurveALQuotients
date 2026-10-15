// Equations produced by the implementation of Molina's algorithm in this directory.
// One representative per curve; each run also returns the other models of the same curve that the
// translation freedom allows, and these are all isomorphic over Q (checked by compare.m).
//
// Rows <D, m, genus, R(x)> for Y = X_0^D(1)/<w_m>  (m = 1: the curve itself), y^2 = R(x).
// The coordinate is NOT Molina's: the points at infinity come from the imaginary quadratic field
// of smallest discriminant that admits them, which is not always his choice, so the models agree
// with his only up to isomorphism.

Px<x> := PolynomialRing(Rationals());

// All ten rows of his Table 2, each checked isomorphic over Q to the published equation.
molina_reproduced := [
  <57,  3, 2, -4*x^6 + 168*x^5 - 2831*x^4 + 24008*x^3 - 103538*x^2 + 190788*x - 52647>,
  <39,  1, 3, -7*x^8 + 102*x^7 - 627*x^6 + 1332*x^5 - 1478*x^4 + 1038*x^3 - 492*x^2 + 144*x - 19>,
  <55,  1, 3, -3*x^8 - 8*x^7 - 40*x^6 - 118*x^5 - 198*x^4 - 184*x^3 - 95*x^2 - 26*x - 3>,
  <35,  5, 2, -8*x^6 + 68*x^5 - 172*x^4 + 8*x^3 + 441*x^2 - 190*x - 375>,
  <51, 17, 2, -7*x^6 - 60*x^5 - 154*x^4 - 84*x^3 + 77*x^2 - 12*x>,
  <65, 13, 2, -7*x^6 + 18*x^5 + 31*x^4 - 46*x^3 - 11*x^2 + 28*x - 8>,
  <65,  5, 2, -7*x^6 + 4*x^5 + 56*x^4 + 2*x^3 - 116*x^2 - 80*x - 15>,
  <69, 23, 2, -4*x^6 - 72*x^5 - 487*x^4 - 1578*x^3 - 2527*x^2 - 1836*x - 504>,
  <85,  5, 2, -3*x^6 + 16*x^5 + 20*x^4 - 50*x^3 - 36*x^2 + 20*x - 7>,
  <85, 85, 2, x^6 - 10*x^5 + 23*x^4 - 10*x^3 - 15*x^2 + 8*x + 4>
];

// Two quotients for which this repository holds no model and which are not in Molina's table.
// Together with X_0^85(1)/<w_5> and X_0^85(1)/<w_85> above, which are in his table but which this
// repository also lacks, that is four of the 43 missing models.
new_equations := [
  <115, 115, 2, x^6 + 2*x^4 - 10*x^3 + 5*x^2 - 6*x + 1>,
  <161, 161, 2, x^6 + 6*x^5 + 17*x^4 + 22*x^3 + 26*x^2 + 12*x + 1>
];
