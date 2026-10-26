# The computation behind `paper/hyperelliptic-21-58.tex`

`X_0^21(58)*` is proved hyperelliptic by counting: four CM divisors of degree two whose differences
are torsion, against a torsion bound of 2. The torsion relations are certified by the dual test
against weight-3/2 forms, in the form the note's lemma gives for a level with `2 || N`: every
Heegner component of the star is `(U(4) + 2)` applied to a Kohnen plus-space form of level 2436
that vanishes wherever the plus-type Heegner series of `X_0^21(29)*` is the zero divisor.

| file | what it does |
|---|---|
| `stepH-plus.m` | computes the cusp forms of weight 3/2 and level 2436 (233 of them; cached under `data/dual/`, not committed), the space `W+` above (dimension 3), `V = (U(4) + 2) W+`, and decides for each target divisor whether its coefficient vector is orthogonal to `V`, which certifies the divisor torsion |

Run from the repository root:

    magma -b DD:=21 NN:=58 MM:=2436 BB:=4000 TGTS:="[[<-16,1>,<-120,1>,<-88,-1>,<-100,-1>]]" tools/proof-21-58/stepH-plus.m < /dev/null

The tests: `tests/Hyperelliptic21_58.m` (a few minutes) checks that the eight points are rational,
identifies the Jacobian with `A_f`, `f = 609.2.a.d`, among the newform orbits Jacquet–Langlands
allows (so that the correspondence to level 29 is injective on it), and recomputes the torsion bound
from the trace formula; `tests/_offline/Hyperelliptic21_58.m` (about half an hour) recomputes the
cusp forms to precision 1600 and runs the test on the four fibre differences and five controls.

Fibres with an odd discriminant (`{-7, -312}`, `{-28, -840}`) are not used: in the pole-at-infinity
dictionary an odd discriminant comes glued to four times itself.
