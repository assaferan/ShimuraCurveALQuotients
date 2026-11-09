# The computation behind `paper/hyperelliptic-77-10.tex`

`X_0^77(10)*` is proved hyperelliptic by counting: five CM divisors of degree two with pairwise
torsion differences, against a torsion bound of 2. The relations are certified by the dual test for
a level with `2 || N` (the lemma of `paper/hyperelliptic-21-58.tex`): every Heegner component of the
star is `(U(4) + 2)` applied to a Kohnen plus-space form of level 1540 vanishing wherever the plus-type
Heegner series of `X_0^77(5)*` vanishes.

| file | what it does |
|---|---|
| `stepH-plus.m` | the cusp forms of weight 3/2 and level 1540 (137 of them, to precision 8000; cached under `data/dual/`, not committed), the space `W+` (dimension 3), `V = (U(4) + 2) W+`, and the decision for each target divisor |
| `scan.m` | the same, run as a search: over every CM divisor of degree at most two on the star in a discriminant `d = 0, 8, 12 mod 16` with `|d| <= 8000` (single rational points, pairs of them, and quadratic points), it lists the ones whose difference from the base fibre is certified torsion. Against `{-16, -100}` it returns `{-280, -660}`, `{-760, -1540}`, `{-840, -1320}` and the quadratic point of discriminant -4180, and nothing else |

Run from the repository root (the scan needs the package attached for the fields of definition):

    magma -b DD:=77 NN:=10 MM:=1540 BB:=8000 TGTS:="[[<-280,1>,<-660,1>,<-16,-1>,<-100,-1>]]" tools/proof-77-10/stepH-plus.m < /dev/null
    magma -b DD:=77 NN:=10 MM:=1540 BB:=8000 CACHE:=data/dual/basis_1540_8000.txt TGTS:="[[<-16,1>,<-100,1>]]" tools/proof-77-10/scan.m < /dev/null

The tests: `tests/Hyperelliptic77_10.m` (a minute) checks the fields of definition of the nine
discriminants, identifies the Jacobian with `A_f`, `f = 385.2.a.e` (level 385, old at 2), and
recomputes the torsion bound; `tests/_offline/Hyperelliptic77_10.m` (about twenty-five minutes)
recomputes the cusp forms and runs the test on the four relations and five controls.
