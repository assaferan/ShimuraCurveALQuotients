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
orders, and at level 924 the products of a rung 75 below the target span it while lower rungs do
not (HANDOFF.md, 2026-10-11).

| level | bound B | shifts | spanning checks (rank of q-expansions at oo, products = enumerated rung) |
|---|---|---|---|
| 204 | 32 (= k) | 1 | (204, 20, 32) x shifts -> (204, 45, 32): 63 = 63 (`tests/TShiftCusp0.m`) |
| 380 | 93 | 58, poles 65..93 | (380, 92, 0) -> (380, 158, 0): 132 = 132 (`tests/TShiftProducts.m`); (380, 28, 65) -> (380, 93, 65): 132 = 132 (campaign `cusp0ladder/ladder.m`) |
| 420 | 104 (= k) | 1, with the reduced rung `tshift_core_420.txt` | (420, 145, 0) -> (420, 250, 0): the ring of forms of X_0^210(1) converges from the 212 products (PROVENANCE.md, level 420) |
| 780 | 274 | 410, poles 197..274 | (780, 274, 0) -> (780, 472, 0): rank 396 = the dimension, from 1049 products kept four per pole order and from 12521 kept 69 per order alike (lovelace, 2026-10-12) |
| 924 | 325 | 282, poles 236..325 | (924, 250, 0) -> (924, 325, 0): 237 = 237, the 565 products against the 8674 enumerated points (lovelace, 2026-10-12) |

The script keeps every product with pole at most the requested `n` while the set has at most
`NMZ_SHIFT_KEEP` points (default 20000), and beyond that as many per pole order as the limit allows,
at every cusp-0 order `m`. A product set that is too small never gives a wrong basis, only one the
search rejects.
