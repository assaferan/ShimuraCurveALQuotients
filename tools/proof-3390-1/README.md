# The computation behind `paper/hyperelliptic-3390-1.tex`

`X_0^3390(1)` modulo `<w_5, w_6, w_113>` is the only genus-four quotient of that curve by an index-two subgroup of
its Atkin-Lehner group, and it is a double cover of the genus-two star branched at the rational CM points of
discriminants -3 and -40. By the normal form of `paper/star-2002-proof.tex` (Lemma 3) it is hyperelliptic exactly when
those two points form a fibre of the star's hyperelliptic map. The Jacobian of the star has no rational torsion, so
the fibre statement is the statement that `P_-3 + P_-40 - P_-43 - P_-67` is torsion, which the dual test of
`paper/hyperelliptic-15-38.tex` (Section 2) decides against the two-dimensional space of cusp forms of weight 3/2 on
`Gamma_0(6780)` vanishing at the indices without CM points.

| file | what it does |
|---|---|
| `branch.m` | the genera, the fixed points of every `w_m` outside `H` by CM discriminant (Ogg), the degrees of the CM places involved, the traces of Frobenius of the star at 7, 11, 13 and the torsion bound of its Jacobian: `magma -b tools/proof-3390-1/branch.m` |
| `basisS.m` | the cusp forms of weight 3/2 on `Gamma_0(M)` to precision `B`, written to `basis_<M>_<B>.txt` (hours at 6780/1100) |
| `withbasis.m` | the dual test of `tests/_offline/Hyperelliptic3390_1.m` run from a saved basis `basis_6780_1100.txt` (the rows are the cusp forms of weight 3/2 on `Gamma_0(6780)` to precision 1100, as written by `basisS.m` in this directory, `magma -b MM:=6780 BB:=1100 tools/proof-3390-1/basisS.m`), for a quick rerun once the basis exists: `magma -b BASIS:=basis_6780_1100.txt tools/proof-3390-1/withbasis.m` |

The tests: `tests/Hyperelliptic3390_1.m` (seconds) checks the genera, that `H` is the unique genus-four subgroup,
the branch locus, the degrees of the CM places, the Hecke data at 7 and the torsion bound;
`tests/_offline/Hyperelliptic3390_1.m` (many hours: it recomputes the basis) certifies the two fibre relations and
refuses two non-fibres.
