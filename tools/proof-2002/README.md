# The computations behind `paper/star-2002-proof.tex`

Each file is one step of the proof that `X_0^2002(1)*` is the stored curve, and each prints the
quantity the proof quotes.  Run from the repository root with `< /dev/null`.

| file | step | what it establishes |
|---|---|---|
| `step0-conductor-euler.m` | 0 | the odd part of the conductor, and the full degree-4 Euler factor at every good prime below 120 against the newform |
| `step0-conductor-at-2.m` | 0 | the exponent 2 of the conductor at 2, by Dokchitser–Doris (needs CHIMP's `Genus2Conductor`; run on lovelace) |
| `step2-endomorphisms.m` | 2 | `End_Q(Jac C) = O_K`, certified (needs CHIMP's `endomorphisms`; run on lovelace); its output is kept as `step2-endomorphisms.log`, the explicit correspondence being the certificate |
| `stepA-orbit.m` | 4 | of the five two-dimensional orbits at level 2002, exactly one has Jac(C)'s Euler factors |
| `stepB-reducibility.m` | 5 | gcd over 24 primes p = 1 (mod 2002) of Norm(a_p − 1 − p) is 1 |
| `stepB-congruence.m` | 5 | gcd of (a_p − σa_p)² is exactly 5 |
| `step5-ramified.m` | 5 | the 5-rank of Jac(C)(F_p) is odd at p = 29, so Jac(C)[5] is not a split extension |

Two of the inputs need packages that are not in this repository: `Genus2Conductor` and
`endomorphisms`, both bundled in `edgarcosta/CHIMP`.  On lovelace that bundle is at
`~/GitHub/CHIMP`; its `MagmaPolred` submodule had to be checked out, and line 41 of
`MagmaPolred/MagmaPolred.m` lacked the semicolon after `end intrinsic`, without which nothing in
`endomorphisms` loads.
