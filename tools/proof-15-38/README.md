# The computation behind `paper/hyperelliptic-15-38.tex`

`X_0^15(38)/<w_2, w_15, w_19>` is a genus-three double cover of the genus-one star `E = X_0^15(38)*`,
branched at the four rational CM points of discriminants -228, -760, -12, -40. It is hyperelliptic
because `P_-228 + P_-760` and `P_-12 + P_-40` are linearly equivalent on `E` exactly, not only modulo
torsion, and because the cover is the untwisted one. The dual test gives the pairing modulo torsion;
the torsion component is read off the component groups of `E` at 3 and 5, where every branch point
reduces to a known component of the Cerednik-Drinfeld special fibre (Molina's theorem on the
specialisation of CM points), and the rational torsion of `E` injects into the component group.

| file | what it does |
|---|---|
| `isogenyclass.m` | the traces of Frobenius of the star from the trace formula at every good prime up to the Sturm bound 240, against 570e4 |
| `specialise.m` | rebuilds the Cerednik-Drinfeld graph at `p = 3` or `5` from one list of ideal classes, counts optimal embeddings of the branch CM orders into the right orders of the vertex and edge classes, forms the quotient by `W` with the thicknesses, and prints the class of `P_-228 + P_-760 - P_-12 - P_-40` (and of the other pairings) in the critical group |
| `realisable.m` | why Borcherds forms cannot do this: against all 113 cusp forms of level 1140 no multiple of the pairing divisor is the divisor of a Borcherds form; the two extra obstruction forms are the Hecke eigen-systems of 570l and 570m |

Run from the repository root (the second needs the branch with `DualGraphMetricQuotientFromData` only
for a cross-check of the thicknesses, which the test does without it):

    magma -b tools/proof-15-38/isogenyclass.m < /dev/null
    magma -b PP:=3 tools/proof-15-38/specialise.m < /dev/null
    magma -b PP:=5 tools/proof-15-38/specialise.m < /dev/null

The tests: `tests/Hyperelliptic15_38.m` (half a minute) identifies the star with 570e4, checks that
its torsion lies off the identity component at 3 and 5, recomputes both component groups, places the
four branch points and shows the pairing divisor has class zero, and compares the dual graph of the
quotient by `H` with the reduction of the two candidate double covers; `tests/_offline/Hyperelliptic15_38.m`
(about six minutes) recomputes the cusp forms of level 1140 and certifies the pairing torsion.
