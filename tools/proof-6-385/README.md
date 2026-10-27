# The computation behind `paper/hyperelliptic-6-385.tex`

`X_0^6(385)` modulo the index-two subgroup of its Atkin-Lehner group is a genus-three double cover of the
genus-one star `E = 2310q1`, branched at the CM places of discriminants -40, -1540 and the place of degree two of discriminant -9240. It is hyperelliptic because the two
rational branch points add up to the place of degree two in the group of `E` exactly, and because the cover is the
untwisted one. The method is that of `paper/hyperelliptic-15-38.tex` and `tools/proof-15-38/`: the dual test
gives the pairing modulo torsion; the component groups of `E` at 2 and 3, where the branch places reduce to known
components of the Cerednik-Drinfeld special fibre (Molina's theorem on the specialisation of CM points), remove the
torsion; and the Prym of the cover, whose L-polynomial at 19 is irreducible over Q, shows the cover is untwisted.

| file | what it does |
|---|---|
| `isogenyclass.m` | the traces of Frobenius of the star from the trace formula at every good prime up to the Sturm bound, against a Cremona curve: `magma -b DD:=6 NN:=385 LAB:=2310q1 tools/proof-6-385/isogenyclass.m` |
| `specialise.m` | the Cerednik-Drinfeld graph at a prime of the discriminant, the quotient by `W` with thicknesses, the critical group, the number of star points of each CM cycle on each component, and the class of a divisor: `magma -b DD:=6 NN:=385 PP:=p DS:="[...]" TGT:="[<d, m>, ...]" tools/proof-6-385/specialise.m` |
| `prym.m` | the L-polynomials of the Prym of the double cover at small primes, factored over Q |

The tests: `tests/Hyperelliptic6_385.m` (a few minutes) identifies the star, checks its torsion against the
identity components, recomputes the component groups, places the branch places, shows the pairing divisor has
class zero and the Prym is simple; `tests/_offline/Hyperelliptic6_385.m` recomputes the cusp forms of weight 3/2
at the level of the dual test and certifies the pairing torsion.
