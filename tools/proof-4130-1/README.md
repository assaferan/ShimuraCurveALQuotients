# The computation behind `paper/hyperelliptic-4130-1.tex`

`Y = X_0^4130(1)*`, the quotient of `X_0^4130(1)` by its full Atkin-Lehner group, has genus three and no further
quotient, so it is hyperelliptic exactly when two distinct divisors of degree two on it are linearly equivalent. Its
Jacobian has no rational torsion (the orders over `F_p` have gcd 1), so it is enough that the difference of two such
divisors be torsion, which the dual test of `paper/hyperelliptic-15-38.tex` (Section 2) certifies for the
differences of the four fibres `P_-35 + P_-532`, `P_-163 + P_-595`, `P_-315 + P_-520`, `P_-340 + P_-420` found by the
Shimizu screen (any one difference suffices).

| file | what it does |
|---|---|
| `data.m` | the genus of the star, the degrees of the eight CM places, the traces of Frobenius at 3, 11, 13 with the characteristic polynomial of `T(p^2)` on `W`, and the torsion bound of the Jacobian: `magma -b tools/proof-4130-1/data.m` |
| `basisS.m` | the cusp forms of weight 3/2 on `Gamma_0(M)` to precision `B`, written to `basis_<M>_<B>.txt` (hours at 8260/1200) |
| `withbasis.m` | the dual test of `tests/_offline/Hyperelliptic4130_1.m` from a saved basis: `magma -b BASIS:=basis_8260_1200.txt tools/proof-4130-1/withbasis.m` |

The tests: `tests/Hyperelliptic4130_1.m` (about four minutes) checks the genus, the rationality of the eight points,
the Hecke data at 3 and the torsion bound; `tests/_offline/Hyperelliptic4130_1.m` (many hours: it recomputes the
basis) certifies the six differences torsion, with a Hecke check on `W` and two refused controls.
