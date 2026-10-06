# The Borcherds obstruction: parity census (2026-10-06)

Borcherds' criterion: a divisor is the divisor of a Borcherds product iff it pairs to zero against
every form of the obstruction space (the weight-3/2 cusp forms of the dual Weil representation). The
search's divisor map is the matrix `coeffs_trunc`, its image the Borcherds divisors reachable from
the weakly holomorphic basis, and the annihilator of that image (`Kernel(Transpose(coeffs_trunc))`,
normalised to primitive integers) is the obstruction space in the CM-discriminant coordinates. A
double cover depends on its branch divisor only mod 2, so a target whose pairings are ALL EVEN can be
shifted by an even divisor E with phi(E) = -phi(target) -- the cover is unchanged and the shifted
divisor is a Borcherds divisor. Odd pairings close that door.

* `obstpair.m`: driver; run from a checkout of branch `obstruction-pairing` (BorcherdsForms.m with
  the `OBSTPAIR=1` print after `found_v`) as
  `OBSTPAIR=1 magma -b DD:=38 NN:=5 .../obstruction/obstpair.m`.
* `parity.py`: tabulates the printed lines per key and anchor triple.
* `pairings_<base>.txt`: the distinct printed lines; `parity_2026-10-06.txt`: the summary.

## Results

| base | failing key | obstruction dim | failing triples | pairing values | parity | gcd(phi) |
|---|---|---|---|---|---|---|
| X_0^38(5) (control) | 11 | 1 | 504 | -2 ... -26 | all even | 1 |
| X_0^14(23) | 12 | 1 | 120 | +-12 | all even | 1 |
| X_0^22(19) | 12 | 2 | 24 | +-8, -12, -14, +-26 | all even (both generators) | 1 |
| X_0^6(109) | 12 | 1 | 60 | 220, +-660 | all even | 1 (phi(3) = 37, phi(4) = -27) |

The control reproduces the 2026-08-29 measurement exactly (key 11, anchor -4 with -19, -11:
pairing -22, phi(4) = 2, phi(760) = -7). X_0^6(109): its plain rerun fails in the search ("Failed to
find all Borcherds forms"), so it is in the obstructed class (the screen's remark "reads 1 0 0 and
builds" is about the deficit reading); the instrumented run took 40 minutes. Its pairings are large,
so the correcting even divisor will have larger coefficients than at the other three.

Each triple appears once per rung of the m-ladder in the printed lines (the while loop re-enters the
triple search), which is why the triple counts are multiples of the number of anchors.

## What the hatch needs in the code

1. In the search: when a target is not in the image, solve the integer system
   sum_d e_d phi_b(d) = -phi_b(target)/2 for every generator b (gcd(phi) = 1 makes the one-generator
   case solvable; the two-generator case at 22_19 is a 2 x n system), add E = sum 2 e_d Z(d) to the
   target, and solve for the form with divisor target + E.
2. Downstream: `assert Set(div_f) eq {...}` after the solve demands the divisor equal the ramification
   divisor exactly, and the Schofer stage builds the cover's equation from the form's values at the CM
   points; both must learn that the form carries extra even components (double zeros and poles at
   the CM points of E), i.e. y^2 = f(t) * (squares). Prefer E >= 0 (extra double zeros only) when the
   system allows it.
3. Validate on X_0^38(5) first (its covers of genus >= 1 against the trace formula), then the three
   genus-0 bases.
