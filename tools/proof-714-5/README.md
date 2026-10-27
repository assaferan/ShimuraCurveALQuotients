# The computation behind `paper/hyperelliptic-714-5.tex`

`X_0^714(5)` modulo the index-two subgroup of its Atkin-Lehner group is a genus-three double cover of the
genus-one star `E = 3570x2`, branched at the CM places of discriminants -340, -840 and the place of degree two of discriminant -14280. It is hyperelliptic because the two
rational branch points add up to the place of degree two in the group of `E` exactly, and because the cover is the
untwisted one. The method is that of `paper/hyperelliptic-15-38.tex` and `tools/proof-15-38/`: the dual test
gives the pairing modulo torsion; the component groups of `E` at 3, 7 and 17, where the branch places reduce to known
components of the Cerednik-Drinfeld special fibre (Molina's theorem on the specialisation of CM points), remove the
torsion; and the Prym of the cover, whose L-polynomial at 11 is irreducible over Q, shows the cover is untwisted.

| file | what it does |
|---|---|
| `isogenyclass.m` | the traces of Frobenius of the star from the trace formula at every good prime up to the Sturm bound, against a Cremona curve: `magma -b DD:=714 NN:=5 LAB:=3570x2 tools/proof-714-5/isogenyclass.m` |
| `specialise.m` | the Cerednik-Drinfeld graph at a prime of the discriminant, the quotient by `W` with thicknesses, the critical group, the number of star points of each CM cycle on each component, and the class of a divisor: `magma -b DD:=714 NN:=5 PP:=p DS:="[...]" TGT:="[<d, m>, ...]" tools/proof-714-5/specialise.m` |
| `prym.m` | the L-polynomials of the Prym of the double cover at small primes, factored over Q |

The tests: `tests/Hyperelliptic714_5.m` (a few minutes) identifies the star, checks its torsion against the
identity components, recomputes the component groups, places the branch places, shows the pairing divisor has
class zero and the Prym is simple; `tests/_offline/Hyperelliptic714_5.m` recomputes the cusp forms of weight 3/2
at the level of the dual test and certifies the pairing torsion.
