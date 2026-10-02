# The Galois-cover route to the obstructed non-squarefree bases

The normaliser of Gamma_0^D(N) exceeds the Atkin-Lehner group exactly when h > 1 (h the largest
divisor of 24 with h^2 | N); see normaliser.m. Each obstructed curve is then a Galois cover of
X_0^D(N/h^2) with group SL_2(F_h), and the star curve is the quotient of that cover by the extra
Atkin-Lehner involution. The star map's ramification at a base star point P is
[G:V] - #orbits(H_P on G/V), computed group-theoretically in branchA4.m.

## Where it works

* `15_4`: degree 3 over X_0^15(1)^*, R = 4. The degree-6 cover is the rigid Belyi map
  lambda -> j, and with branch values read off models_15_1.m it reproduces Tu's t_4 to 60
  digits (cover15.m).
* `21_4`: degree 3 over X_0^21(1)^*, R = 4, four simple branch points at the star points of
  d = -4, -28 and the two of -84. A UNIQUE degree-3 map has those critical values, it is defined
  over Q although two of them are not, and its fibres carry the right CM fields
  (hurwitz21.m, fibres21.m):  1/s = (t^3 - 4/3 t + 16/27)/(t^2 + 29/12 t + 22/9).
* `10_9`: degree 6 over X_0^10(1)^*, R = 10, and the base's own special points account for it
  exactly: the d = -3 elliptic point has stabiliser image C_3 (contributing 4) and the three
  w-fixed points d = -40, -20, -8 have image C_2 (2 each).

## Where the enumeration is INCOMPLETE -- measured, not suspected

Riemann-Hurwitz does NOT close for three of the bases using the base's elliptic points and
w-fixed points alone:

    14_9   degree 6, needs R = 10; two d = -56 points, d = -8 and the d = -4 elliptic point give
           at most 2+2+2+3 = 9, and no subgroup of A_4 contributes the missing 1 -- so at least
           one branch point is missing (one more C_2 point would give exactly 10).
    22_9   degree 6, needs R = 12; -88, -11, -4, -3 give 2+2+2+4 = 10; two short.
    15_8   degree 4 over X_0^15(2)^*, needs R = 6; its five base points contribute at most 1
           each under C_2 x C_2, so at most 5; at least one short.

The missing points are the fixed points of the EXTRA normaliser elements themselves. On a Shimura
curve there are no cusps, so the element that is parabolic in the classical picture
((1 1/h; 0 1) for Gamma_0(h^2 N')) must be elliptic here, and its fixed points are CM points of
the discriminant t^2 - 4n(alpha) that do not appear among the base's w-fixed points. Computing
those discriminants is the next step; until it is done, the branch data of 14_9, 22_9 and 15_8 is
not known, and the three cases above are not solvable by this route. ⇒ 15_4, 21_4 and 10_9 close
on their own, which is evidence for them but not a proof that nothing was missed.
