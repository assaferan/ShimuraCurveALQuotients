# Candidate models of the genus-two Atkin–Lehner stars, from the periods

⚠⚠ **READ THIS FIRST. THESE ARE CANDIDATES, NOT ESTABLISHED MODELS.** The point counts do not
identify the curve. All the principal polarizations on a modular abelian surface give curves whose
Jacobians are isomorphic as *unpolarized* varieties, so they share an L-function and therefore
agree in #C(F_p) at **every** prime. A count-based check confirms the isogeny class and is blind to
which curve one has. See the demonstration below.

## The demonstration, at level 35

X_0^35(1)/w_5 has an independent model from the Guo–Yang route, stored in
`data/models/models_35_1.m`. The period route returns a different curve:

    period route : y^2 + (x^3 + x) y = -x^5 - x^4 + 3x^3 - 3x^2 - 5x
    Guo-Yang     : y^2 = -7/256 x^5 - 25/64 x^4 - 277/128 x^3 - 345/64 x^2 - 7/256 x

    G2 invariants, period route : [ -21507498573824/52521875, -342847684608/10504375,
                                    -5511792896/2100875 ]
    G2 invariants, Guo-Yang     : [ 187369943065427968/102942875, 396126074044416/2100875,
                                    -27001865216/1715 ]

The invariants differ, so the two are not isomorphic even over the algebraic closure — and their
point counts agree at every prime tested (3, 11, 13, 17, 19, 23, 29, 31, all eight). This is the
phenomenon of González–Guàrdia–Rotger: non-isomorphic genus-two curves sharing an unpolarized
modular Jacobian. It is also a diagnosis: at level 35 the search over small combinations of the
two Néron–Severi generators found only one principal polarization class, and the Guo–Yang curve
belongs to another that the search did not reach.

## The candidates

Produced by `periodstar.m` at 150 digits. Each has the right Jacobian up to isogeny; whether it is
the star is **not** established.

**X_0^2002(1)\***

    y^2 + (x^2 + x) y = 16x^6 + 4x^5 + 180x^4 - 456x^3 + 196x^2 + 64x + 4

Cardona–Quer invariants recognised with margins 2.9e92, 1.8e101, 4.0e109:

    4443541864293086393626096432993/8718043021712850681856
    -34439169873536221631669217/2177333422006206464
    186506522749409116433/543789565935616

**X_0^2030(1)\***

    y^2 + (x^2 + x) y = x^6 + 29x^5 + 13x^4 + 9x^3 + 159x^2 + 149x + 39

    1200089628966417104818086601/3102846490240000000000
    123175783384857423892659/8865275686400000000
    63643057903538421401/25329359104000000

Each matches p + 1 − a_p untwisted at all 13 primes below 60 of good reduction, with discriminant
supported exactly on the primes of D — {2,7,11,13} and {2,5,7,29}. Both facts are true of the star
and were not imposed, so they are evidence that the isogeny class is right; neither distinguishes
the curve from the others sharing that Jacobian.

## What a real selector would have to do

Distinguish principal polarizations, which point counts cannot. Two candidates:

* **The CM points.** The star carries CM points whose fields of definition are known
  arithmetically; on the correct curve the rational ones are rational. In the analytic coordinate
  x = f2/f1 of the Shimura side, rational points lie on a common circle or line, which is a
  computable test. It is validated at X_0^35(1)/w_5, where the CM points of discriminants −28, −43
  and −67 (all of class number one) come out collinear to five digits, and it fails at the 2002
  candidate — consistent with that candidate not being the star.
* **Enumerating the polarizations properly.** At level 35 the correct curve is evidently missed by
  the search over coefficients in [−6, 6]; a complete enumeration of the principal polarizations
  up to isomorphism, in the sense of González–Guàrdia–Rotger, would produce the full candidate list
  and the Guo–Yang curve should appear in it.

## Status

Level 23 reproduces `SmallModularCurve(23)` exactly and level 67 returns a curve that all of its
polarizations agree on, so the route is not generally wrong; what it lacks is any means of knowing
when it is. Until a selector exists, nothing here should be used as a model of a star.
