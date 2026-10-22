# The polytope cache and the shift sets

`polymake_solution_M_n_m` holds the lattice points of the polytope `(M, n, m)`: the exponent
vectors of the weight-1/2 eta quotients of level `M` with a pole of order at most `n` at the cusp
oo and at most `m` at the cusp 0 (the "rung"). The Borcherds search reads a cached file as is, and
otherwise asks `nmzsolve.py` for it (`polymake_script_*` are the inputs of the retired polymake
backend, kept as a record of which rungs were requested).

A file whose first line is a comment `// derived: ...` is not the lattice points but a spanning
set of the rung's forms, built by the solver script as the products of a lower enumerated rung with
the level's weight-0 shift set (below). The search checks the rank of what it receives, and
enumerates a rung whose derived set falls short.

## The shift sets `tshift_w0_M.txt`

The weight-0 eta quotients of level `M` with poles at oo only and pole order at most `B`,
written by the solver script itself, `python3 nmzsolve.py M B 0 polymake/tshift_w0_M.txt 0 1 0`.
Multiplying a rung `(M, n, m)` by them gives forms of the rung `(M, n + B, m)` with the same pole
at 0, so a tall rung is reached from a lower one in `ceil((n_tall - n_low) / B)` layers of
products, with no Normaliz run. The smallest pole order of a weight-0 quotient at the level is the
pole `k` of the Hauptmodul-like quotient the search uses as its climbing step, so any bound
`B >= k` reaches every rung; the bound used at each level is the pole order of the first rung the
search asks for there (or `k` itself), and whether the products span the tall rung is a measured
fact, checked by the rank of the q-expansions, not a theorem: eta quotients are sparse at low pole
orders, and at level 924 the products of the rung (924, 250, 0) span (924, 325, 0) while the
products of the rungs at 120, 160, 200 and 230 do not (ranks of the q-expansions of each product
set against the enumerated rung, 2026-10-11).

The spanning checks in the table compare the rank of the q-expansions at oo of a set of exponent
vectors, computed by `tools/etarank.m` over GF(p) for two large primes, independently of the
repository's code, with the **dimension** of the space the rung is meant to span. For eta quotients
of level M = 4*M_0 with pole order at most n at the cusp oo and none at the other cusps, that
dimension is

    n + sum over d | M_0 of floor(d/4) + 1 - g(X_0(M)),

which is 396 at `(780, 472, 0)` and 237 at `(924, 325, 0)`: g(X_0(780)) = 157 and g(X_0(924)) = 181,
both of index 2016 and 2304 with 24 cusps and no elliptic points. A rank over a finite field is at
most the rank over Q, which is at most the dimension, so a measured rank *equal to the dimension* is
conclusive, where equality of two measured ranks would only be two measurements agreeing.

At 780 the set is the one the solver itself serves for the rung `(780, 472, 0)`, formed from the
committed rung 274 and the 410 shifts with no Normaliz run:

    python3 nmzsolve.py 780 472 0 /tmp/sp780
    magma -b M:=780 file:=/tmp/sp780 npole:=472 tools/etarank.m < /dev/null

and its rank is 396 at both primes, the dimension. Those 12521 vectors are a thinned selection and
not every product: the shifts give 333051 points of pole order at most 472, of which the script
keeps 69 per pole order under the `NMZ_SHIFT_KEEP` cap described below, fewer even than the 12729
points of the base rung itself. At 924 the enumerated rung `(924, 325, 0)` is committed, and

    magma -b M:=924 file:=polymake/polymake_solution_924_325_0 npole:=325 tools/etarank.m < /dev/null

gives 237 on its 8674 vectors, again the dimension; the 565 products of the rung `(924, 250, 0)`
with the shifts have rank 237 as well, but they were formed on a machine whose copy of that lower
rung is not in the repository, so that half of the row is a measurement rather than a reproducible
check. The 204 and 380 checks are the two tests named in the table.

| level | bound B | shifts | spanning checks (rank of the q-expansions at oo against the dimension) |
|---|---|---|---|
| 204 | 32 (= k) | 1 | (204, 20, 32) x shifts -> (204, 45, 32): 63 = 63 (`tests/TShiftCusp0.m`) |
| 380 | 93 | 58, poles 65..93 | (380, 92, 0) -> (380, 158, 0): 132 = 132 (`tests/TShiftProducts.m`); (380, 28, 65) -> (380, 93, 65): 132 = 132 (`vvdata/weyl-campaign/cusp0ladder/ladder.m` on the campaign branch) |
| 420 | 104 (= k) | 1, with the reduced rung `tshift_core_420.txt` | (420, 145, 0) -> (420, 250, 0): the ring of forms of X_0^210(1) converges from the 212 products (`data/models/PROVENANCE.md`, the level 420 section) |
| 780 | 274 | 410, poles 197..274 | (780, 274, 0) -> (780, 472, 0): the 12521 vectors the solver serves, a thinned selection of 69 per pole order out of 333051 products, have rank 396 = the dimension (`tools/etarank.m`, commands above) |
| 924 | 325 | 282, poles 236..325 | (924, 250, 0) -> (924, 325, 0): the committed enumerated rung has rank 237 = the dimension (`tools/etarank.m`); the 565 products had rank 237 too, from a lower rung not committed |

The script keeps every product with pole at most the requested `n` while the set has at most
`NMZ_SHIFT_KEEP` points (default 20000), and beyond that as many per pole order as the limit allows,
at every cusp-0 order `m`. A product set that is too small never gives a wrong basis, only one the
search rejects.
