# Saved weight-3/2 bases for the dual test

The dual test (`DualTest.m`) and the certificates of `data/certificates/` read bases of the space of cusp forms
of weight 3/2 on Gamma_0(M) with the theta^3 multiplier, M = 4 * oddpart(D N). They are too large for git (up
to 330 MB each, about 1 GB in all) and are deposited on Zenodo; `bases.md5` in this directory lists every file
with its MD5 checksum and size, and is what the code checks against.

**Format.** `basis_<M>_<B>.txt` is a Magma matrix literal (`Matrix(Rationals(), r, B, [...])`), one cusp form per
row, the coefficients a(1), ..., a(B) of its q-expansion in the columns, r = dim S_{3/2}(Gamma_0(M)). Read it
with `ReadWeightThreeHalvesBasis(M, B)`, which also checks the checksum when one is given.

**How they were computed.** `CuspidalSubspace(HalfIntegralWeightForms(Gamma0(M), 3/2))` and `Basis(S, B + 2)` in
Magma 2.29 on lovelace (MIT), single-threaded, 2026-10-26 to 2026-11-08: minutes below level 2000, hours to a
day at levels 3000-5500, two to four days at levels 7140-8260 (100-125 GB of memory). The same call recomputes
any of them (`ReadWeightThreeHalvesBasis(M, B : Compute := true)`, or `tools/proof-4130-1/basisS.m`).

**Where the code looks.** In `$DUALW_CACHE`, default `~/dualw`. `tools/fetch_dual_bases.sh` downloads the
deposit there and checks every file against `bases.md5`.

**Which certificate uses which basis.** Each certificate names its basis as `<M, B, variant, md5>` in
`c["dualtests"]`; `tests/_offline/DualTestReproduces.m` reproduces the dual tests of the original proofs from them.
