# Heights of the branch points on an elliptic quotient

Scripts run from a checkout of `main` carrying `ExactDoubleCover.m` (PR #123), 2026-10-18.

* `h2.m`: X_0^35(1) and X_0^51(1). The canonical heights of the branch points of X -> E = X/w on
  E are exactly rational multiples of L(E,1) L'(E^d,1) sqrt|d| / Area(E) (Area = covolume of the
  period lattice): 2 and 1/2 at D = 35 (d = -7; the -28 point is -2 times the -7 point, a_2 = 0),
  3/2 and 0 at D = 51 (d = -3; one point is torsion). Thirty digits.
* `scan.m`, `scan.log`: the same over every stored genus-3 hyperelliptic quotient whose completed
  octic is even and whose elliptic quotient has a small rational point: 70 cases. The ratio is
  rational in every case where it is defined (L'(E^d,1) ≠ 0), ranging from 1/4320 to 18, with the
  square of the unit index visible at d = -3 (9/32, 9/14, 9/22, 9/26, 18 = 2·9). Branch points over
  REAL quadratic fields and rational branch points occur (then the comparison quantity must be
  built from L'(E,1) and the twist differently; not done). The rational factor is not yet derived:
  it should be the Gross-Zagier/Zhang constant for the parametrisation X_0^D(N) -> X/W -> E, with
  the unit index, the degree 2|W|, the conductor of the CM order, and the quaternionic Petersson
  norm. Deriving it and checking it on these 70 cases is what turns the heights into an a priori
  bound, which is what the exact route to the level-5460 cases (210_13) needs.
