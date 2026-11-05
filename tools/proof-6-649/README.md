# The computation behind `paper/hyperelliptic-6-649.tex`

`Y = X_0^6(649)*`, the quotient of `X_0^6(649)` by its full Atkin-Lehner group, has genus three and no further
quotient, so it is hyperelliptic exactly when two distinct divisors of degree two on it are linearly equivalent. The
three candidates `P_-24 + P_-132`, `P_-40 + P_-187`, `P_-52 + P_-660` (fibres of the hyperelliptic map found by the
Shimizu screen, which the proof does not use) have pairwise differences that the dual test of
`paper/hyperelliptic-15-38.tex` (Section 2) certifies torsion. The rational torsion of the Jacobian has order dividing
6; its 3-part injects into the component group at 3 (totally degenerate reduction, a torsion-free formal group over
`Z_3`, and a split torus whose `F_3`-points have no 3-torsion), where all three differences have class zero by
Molina's description of the reduction of CM points; so the differences lie in the 2-part, of order at most two, and
two of the three divisors coincide in the Jacobian.

| file | what it does |
|---|---|
| `specialise.m` | the Cerednik-Drinfeld graph at a prime of the discriminant, the quotient by `W` with thicknesses, the critical group, the number of star points of each CM cycle on each component, and the class of a divisor: `magma -b DD:=6 NN:=649 PP:=3 DS:="[-24,-132,-40,-187,-52,-660]" TGT:="[<-24,1>,<-132,1>,<-40,-1>,<-187,-1>]" tools/proof-6-649/specialise.m` (at 3: `Z/102`; at 2: `Z/3 x Z/3 x Z/90`; every difference of the three fibres has class zero at both) |
| `basisS.m` | the cusp forms of weight 3/2 on `Gamma_0(M)` to precision `B`, written to `basis_<M>_<B>.txt` (hours at 7788/1200) |
| `withbasis.m` | the dual test of `tests/_offline/Hyperelliptic6_649.m` from a saved basis: `magma -b BASIS:=basis_7788_1200.txt tools/proof-6-649/withbasis.m` |

The tests: `tests/Hyperelliptic6_649.m` (about ten minutes) checks the torsion bound, the graphs and component
groups at 3 and 2, that the reduction is totally degenerate, and the classes of the three differences;
`tests/_offline/Hyperelliptic6_649.m` (many hours: it recomputes the basis) certifies the three differences torsion,
with a Hecke check on `W` and two refused controls.
