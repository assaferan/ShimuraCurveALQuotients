# The computation behind `paper/hyperelliptic-6-323.tex`

`X_0^6(323)` modulo the index-two subgroup `<w_6, w_19, w_34>` of its Atkin-Lehner group is a genus-three double
cover of the genus-one star `E = 1938h1`, branched at a single place of degree four: the CM points of discriminant
-7752. The four points are permuted cyclically by the class group of `Q(sqrt -1938)` modulo its 2-torsion, the only
Galois-stable pairing is `{P, s^2 P} | {s P, s^3 P}`, and the curve is hyperelliptic because the two halves are
linearly equivalent on `E` and the cover is the untwisted one. The pairing divisor is a CM cycle twisted by the genus
character `(-8, 969)`: it is torsion by the Gross-Zagier formula of Yuan-Zhang-Zhang (the twisted L-function vanishes
to order at least three), and zero by reduction at 2, where the four points reduce two to each of two components whose
difference has order 14 in the component group, so that no distribution of the pairs over the components other than
the ones giving zero is compatible with a torsion divisor. The Prym of the cover has an irreducible L-polynomial at
7, so the cover is untwisted. The method otherwise follows `tools/proof-15-38/` and `tools/proof-6-385/`.

| file | what it does |
|---|---|
| `isogenyclass.m` | the traces of Frobenius of the star from the trace formula at every good prime up to the Sturm bound, against a Cremona curve: `magma -b DD:=6 NN:=323 LAB:=1938h1 tools/proof-6-323/isogenyclass.m` (1938h1 and 1938h2 agree; the component groups at 2 and 3 tell them apart) |
| `twists.m` | the curves of the class 1938h with their torsion over Q and over Q(sqrt -2), their component groups, whether the 2-torsion point lies on the identity component at 2 and 3, and the root numbers and ranks of the twists of E by -8 and 969 |
| `specialise.m` | the Cerednik-Drinfeld graph at a prime of the discriminant, the quotient by `W` with thicknesses, the critical group, and the number of star points of a CM cycle on each component: `magma -b DD:=6 NN:=323 PP:=2 DS:="[-7752]" TGT:="[<-7752,0>]" tools/proof-6-323/specialise.m` (at 2: two points on each of two components at distance 3 on the cycle of length 14; at 3: the critical group has order 3) |
| `prym.m` | the L-polynomials of the Prym of the double cover at small primes, factored over Q |

The test `tests/Hyperelliptic6_323.m` (about two minutes) identifies the star, checks the class group and the genus
character, the root numbers and ranks behind the torsion statement, the torsion of the minus part over both
candidate quadratic fields, the position of the 2-torsion point, the specialisation of the branch place at 2 in the
component groups over Q_2 and over Q_2(sqrt -2), and the irreducibility of the Prym's L-polynomial at 7.
