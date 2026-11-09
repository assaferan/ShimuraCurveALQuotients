# The computation behind `paper/hyperelliptic-39-34.tex`

`X_0^39(34)*` is proved hyperelliptic by counting: five CM divisors of degree two with pairwise
torsion differences, against a torsion bound of 4. The relations are certified by the dual test for a
level with `2 || N` (the lemma of `paper/hyperelliptic-21-58.tex`), in two forms: the pole-at-infinity
dictionary of the star against `(U(4) + 2) W+`, and, for the fibre in the odd discriminant -15, the
plus-type dictionary of `X_0^39(17)*` against `W+` after pushing the divisor down by the
correspondence of the lemma.

| file | what it does |
|---|---|
| `stepH-plus.m` | the cusp forms of weight 3/2 and level 2652 (245 of them, to precision 4000; cached under `data/dual/`, not committed), the space `W+` (dimension 3), `V = (U(4) + 2) W+`, and the decision for each target divisor in even discriminants |
| `stepH-odd.m` | the same test for divisors in odd discriminants: pushes the divisor to `X_0^39(17)*` with the fibre count of the lemma, writes it in the plus-type dictionary, and tests against `W+` |
| `scan.m` | `stepH-plus.m` run as a search over every CM divisor of degree at most two of the star in a discriminant `d = 0, 8, 12 mod 16`; against `{-84, -132}` it returns the quadratic points of discriminants -3432 and -1560 and the pair `{-1092, -372}` |

Run from the repository root:

    magma -b DD:=39 NN:=34 MM:=2652 BB:=4000 TGTS:="[[<-1092,1>,<-372,1>,<-84,-1>,<-132,-1>]]" tools/proof-39-34/stepH-plus.m < /dev/null
    magma -b DD:=39 NN:=34 MM:=2652 BB:=4000 CACHE:=data/dual/basis_2652_4000.txt TGTS:="[[<-15,1>,<-52,1>,<-84,-1>,<-132,-1>]]" tools/proof-39-34/stepH-odd.m < /dev/null
    magma -b DD:=39 NN:=34 MM:=2652 BB:=4000 CACHE:=data/dual/basis_2652_4000.txt TGTS:="[[<-84,1>,<-132,1>]]" tools/proof-39-34/scan.m < /dev/null

The tests: `tests/Hyperelliptic39_34.m` (a minute) checks the fields of definition of the eight
discriminants, identifies the Jacobian with `A_f` of level 663 (old at 2), and recomputes the torsion
bound; `tests/_offline/Hyperelliptic39_34.m` (about thirty-five minutes) recomputes the cusp forms and
runs both forms of the test on the four relations and five controls.
