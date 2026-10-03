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

## The missing points identified, 2026-10-02 (extrapts.m)

The extra normaliser element is the translation by 1/h, the matrix (h 1; 0 h) of reduced norm h^2
and trace t. Classically t = 2h and it is parabolic, fixing a cusp; a Shimura curve with D > 1 is
compact, so the global element is ELLIPTIC and its fixed points are CM points of discriminant
t^2 - 4h^2 with |t| < 2h. Those discriminants are not of the form -q or -4q for q | D N', so such a
point is fixed by no Atkin-Lehner involution and was absent from the enumeration.

For an image of order 2 in the Galois group the element squares into Q^* O^*, forcing t = 0 and
d = -4h^2: **d = -36 for h = 3, d = -16 for h = 2.**

### It closes the three failures, and predicts the one that needed nothing

    14_9   degree 6, R = 10:  -56 (x2), -8, -4 give 2+2+2+2 = 8, and the base has ONE star point
                              of d = -36 -> +2 = 10  EXACTLY
    22_9   degree 6, R = 12:  -88, -11, -4, -3 give 2+2+2+4 = 10, one d = -36 star point
                              -> +2 = 12  EXACTLY
    10_9   degree 6, R = 10:  closed already -- and the base has NO d = -36 points, so the
                              hypothesis demands no extra branch point there.  CONSISTENT.
    15_8   degree 4, R = 6:   its five base points contribute at most 1 each, so it needs a sixth
                              branch point; d = -16 has no points on 15_2, and the only candidate
                              among t^2 - 16 that does is d = -7 (t = +-3), one star point.

### ⚠ What is still open: the trace t at h = 2

The three h = 3 cases are settled by t = 0. The h = 2 cases are not, and the two choices disagree:
15_4's cover is branched at three points only (cover15.m reproduces Tu's values there, so its data
is complete), which needs the extra element to have NO fixed point on 15_1 -- true for d = -16,
false for d = -7, since 15_1 carries a -7 star point. 15_8 needs the opposite. So the actual trace
has to be COMPUTED, per base, not guessed: enumerate the elements of reduced norm h^2 in the Eichler
order, keep those conjugating the order to itself and not of the form (scalar) x (unit), and read
off t. Until that is done, 15_8 and 33_4 are not solvable by this route, while 14_9 and 22_9 are
(their extra point is the d = -36 one).

⚠ **That enumeration is not a short-vector search.** B is INDEFINITE here (split at infinity, as a
Shimura curve needs), so the reduced-norm form is indefinite: the elements of norm h^2 are infinite
in number and `ShortVectors` does not apply. The two routes are (a) the local embedding theory at
the prime dividing h -- which orders of discriminant t^2 - 4h^2 embed optimally in the Eichler order
of level h^2 N' -- or (b) Magma's two-sided ideal machinery, and ⚠ `TwoSidedIdealClassGroup` FAILS
its positive control here (it returns 1 at 15_4, where the normaliser is demonstrably larger), so
its output must not be quoted. Note also that the extra element's fixed points need NOT create a new
branch point: their discriminant can coincide with one that already branches, which is one way 15_4
could need nothing new, and distinguishing that from "no fixed points at all" is part of what the
computation has to settle.

## ⚠⚠ RETRACTION of the group-theoretic framework above (topal.m, 2026-10-02)

**The Galois group is not SL_2(F_h).** Atkin-Lehner-Newman gives the normaliser of Gamma_0(N)
modulo Gamma_0(N) and the Atkin-Lehner involutions as CYCLIC of order h, generated by
tau -> tau + 1/h, not SL_2(F_h) of order 6 or 12. What is SL_2(F_h) is the MONODROMY group of the
degeneracy map X_0^D(h^2 N') -> X_0^D(N') -- S_3 for h = 2, acting on 3 sheets -- which is a
different statement, and the coset counting in branchA4.m used the wrong group for h = 3.
⇒ **10_9's apparent closure is not evidence**, and neither is 33_4's; the contributions table there
(|H| = 2 gives 2, |H| = 3 gives 4, ...) applies to a group that is not acting.

**And a second omission, which is what topal.m measures.** The Hall divisors of D*N include ones
that are not Hall divisors of D*N' -- 4, 12, 20, 60 when N = 4; 9, 18, 45, 90 when N = 9; 8, 24, 40,
120 when N = 8 -- so the TOP curve has Atkin-Lehner involutions the base does not, and their fixed
points (discriminants -m and -4m) are branch points of the map. Those discriminants are fixed by no
involution of the base, which is exactly why the enumeration missed them. But taking all of them
OVERSHOOTS: at 14_9 the new discriminants -36 (w_9), -72 (w_18) and -504 (w_126, two star points)
would add 10 to a required total of 10, on top of the 6 the base's own points already give. So some
of these points are not branch points, or not with the multiplicity the wrong group predicted, and
the earlier "d = -36 closes 14_9 and 22_9 exactly" was arithmetic in the wrong framework.

### What survives

* The normaliser theorem (`normaliser.m`): classical Atkin-Lehner-Newman plus strong approximation.
  Independent of all of this.
* `15_4`: the cover is the Belyi map lambda -> j and `cover15.m` reproduces Tu's six values to 60
  digits from `models_15_1.m`. Verified numerically, not inferred from the framework.
* `21_4`: the degree-3 map 1/s = (t^3 - 4/3 t + 16/27)/(t^2 + 29/12 t + 22/9). Its branch values
  came from the discredited framework, but the solve's OUTPUT carries three checks that do not:
  the map is unique, it is defined over Q although two critical values are not, and its fibres
  generate exactly Q(sqrt(-7)) over the d = -7 point and Q(sqrt(-3)) over the two d = -84 points.
  So the result is probably right and the derivation of its input needs redoing.

### What the branch data actually needs

The local ramification of the degeneracy map X_0^D(h^2 N') -> X_0^D(N') at each elliptic and CM
point, combined with the two W-quotients -- i.e. the classical ramification of X_0(h^2 N') -> X_0(N')
transported to the quaternionic setting, where there are no cusps. That is a bigger computation than
coset counting and has not been done. Until it is, only 15_4 and 21_4 have branch data worth
trusting, and the other five bases are open.

## The degeneracy map's ramification, VERIFIED (degen.m, 2026-10-02)

Rebuilding the branch data from the bottom, the first step is the map before any Atkin-Lehner
quotient. H/Gamma_0(N) -> H/Gamma_0(N') ramifies exactly where the point stabiliser shrinks, i.e.
over the base's ELLIPTIC points, and an elliptic point of order e downstairs has n_top/n_bot
elliptic preimages (index 1) with the rest in orbits of size e, so

    R = sum_{e in {2,3}} n_bot(e) * (deg - n_top(e)/n_bot(e)) * (e-1)/e,   deg = psi(N)/psi(N').

**Riemann-Hurwitz reproduces the genus formula's g_top in all eleven cases tested**: the seven
obstructed bases, the two unobstructed ones, and two squarefree controls (`10_3`, `14_3`).

⚠ The first version of this law assumed p^2 | N kills the elliptic points upstairs. That is true
only for p = 2 (order 2) and p = 3 (order 3) -- exactly the h | 24 condition -- and it FAILED at
`6_25` and `6_49`, where 5 and 7 are split or inert in Q(i) and Q(sqrt -3) and the elliptic points
persist or multiply (`6_25`: 2 of order 2 downstairs, 4 upstairs). With the n_top term restored all
eleven close.

### What is still missing for the star map

Passing from X_0^D(N) -> X_0^D(N') to X^*(D,N) -> X^*(D,N') needs the two Atkin-Lehner quotients,
and that is the step this file retracted above: the containment of the two groups, and how a W-orbit
upstairs sits over a W-orbit downstairs at a point with nontrivial stabiliser. The degeneracy law
above is independent of it and settled; the star step is not, so the branch data of 15_8, 10_9,
14_9, 22_9 and 33_4 remains open.

## ⚠⚠ RETRACTION of the 21_4 hauptmodul, and the structural reason the route stops at 15_4
## (whichobject.m, 2026-10-02)

### ⚠⚠ THE RETRACTION BELOW IS ITSELF WRONG (2026-10-03, Sachi's review of PR #67)

The "no map" statement below is about the projection z -> z, which indeed does not descend. The
DEGENERACY map z -> q z does: conjugating by diag(q,1) carries w_{q^2} into q times an element of
the lower-level group, exactly as X_0(4)/w_4 -> X(1). Checked classically at level 84 = 4*21:
diag(2,1) * (4 -1; 84 -20) * diag(2,1)^-1 = (4 -2; 42 -20) = 2 * (2 -1; 21 -10), and (2 -1; 21 -10)
is in Gamma_0(21) (scratch, 2026-10-03). So X^*(21,4) -> X^*(21,1) EXISTS, of degree
psi(4)/psi(1) * #W'/#W = 3, and "the degree 3 was an arithmetic ratio, not the degree of a map" is
false: it is the degree of this map. ⇒ **The degree-3 Hurwitz solve of hurwitz21.m was about a
real object, and the hauptmodul 1/s = (t^3 - 4/3 t + 16/27)/(t^2 + 29/12 t + 22/9) is NOT retracted;
it is TO BE RECHECKED** -- in particular whether the branch data fed to the solve is the branch data
of this map, and whether its CM fibres (Q(sqrt -7), Q(sqrt -3)) are those of X^*(21,4). The point
about Galois-stable input not being independent confirmation still stands, so the recheck needs an
outside value (a CM value or a point count), not the solve's own rationality. The quotient map by
U = Hall(DN) cap Hall(DN') described below is a different map and everything said about it stands.

**There is no map X^*(D,N) -> X^*(D,N') [WRONG for the degeneracy map; see above].** Such a map needs the top Atkin-Lehner group inside
<Gamma_0(N'), W'>, and the new involutions are not: w_4 has reduced norm 4, so it would have to be
2 times a norm-one unit of the maximal order, i.e. w_4/2 integral, which a primitive element of norm
4 is not. Classically w_4 = (0 -1; 4 0)/2 is visibly outside SL_2(Z), and the reason is geometric --
w_4 sends (E, C) to (E/C, ...), so the forgetful map does not descend. The degree 3 I computed for
"X^*(21,4) -> X^*(21,1)" was an arithmetic ratio psi(N)/psi(N') * |W'|/|W_top|, not the degree of a
map. ⇒ **The degree-3 Hurwitz solve of hurwitz21.m is about an object that does not exist, and the
hauptmodul 1/s = (t^3 - 4/3 t + 16/27)/(t^2 + 29/12 t + 22/9) is RETRACTED.** Its unique, rational
solution with fibres over Q(sqrt -7) and Q(sqrt -3) was a solution of that Hurwitz problem; the
rationality follows from the input being Galois-stable, so it was not the independent confirmation I
took it for.

**What does map down** is X_0^D(N)/U -> X_0^D(N')/U with U = Hall(DN) cap Hall(DN') -- the largest
Atkin-Lehner group acting on both levels -- of degree psi(N)/psi(N'). Tu's t_4 is exactly this for
15_4: he quotients by W_15 = {1,3,5,15}, not by Hall(60). And the genus of that curve decides
whether a hauptmodul exists at all:

    base   U = Hall(DN) cap Hall(DN')   degree   g(X_0^D(N)/U)
    15_4   {1,3,5,15}                  6        0     <- a hauptmodul exists
    21_4   {1,3,7,21}                  6        1
    33_4   {1,3,11,33}                 6        3
    15_8   {1,3,5,15}                  4        1
    10_9   {1,2,5,10}                 12        1
    14_9   {1,2,7,14}                 12        1
    22_9   {1,2,11,22}                12        3

⇒ **15_4 is the ONLY one of the seven whose curve has genus 0.** That is a complete structural
explanation of why Tu treats 15_4 and no other non-squarefree case, and why Guo-Yang's Remark 39
cites him for exactly this curve. **The 15_4 route does not generalise**, and the reason is not
missing branch data or a missing Hauptmodul in the literature: for the other six the object is a
curve of genus 1 or 3, so there is no Hauptmodul to find. Any equation for them has to come from a
positive-genus model, i.e. from the pipeline's own machinery once it handles non-squarefree levels,
not from a Belyi-type solve.

### Everything in this file that still stands

* The normaliser theorem (`normaliser.m`) -- independent of all of the above.
* The local lattice structure at p^e and the m = 0 analysis at a p^2-scaled plane, with its counting
  check (`../level-p2/`) -- independent, and it is what `6_25` and `6_49` need.
* The degeneracy map's ramification law (`degen.m`), verified against the genus formula on eleven
  cases.
* `15_4`: `cover15.m` reproducing Tu's six values to 60 digits, which is about the right object
  (degree 6 out of X_0^15(4)/W_15, genus 0).
