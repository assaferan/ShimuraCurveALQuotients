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

The genus-four quotient `Y = X_0^2002(1)/<w_2, w_13, w_77>` (Section 3 of the note) is a double
cover of the star, and three more scripts carry the star's proof over to it:

| file | what it establishes |
|---|---|
| `stepH-branch.m` | the fixed points of every `w_m` with `m` outside `H`, by CM discriminant: eight of discriminant −11, eight of −91, none else; with Riemann–Hurwitz, the branch locus is one rational point of each |
| `stepH-enumerate.m` | over every fibre of the star above a rational point of height ≤ 10^6 and every twist supported on {−1, 2, 7, 11, 13}, which double covers have the point counts of `Y`: exactly one, the fibre at infinity with `c = 2` |
| `stepH-verify.m` | the resulting curve `Y_2`, its two genus-two quotients (the star and `C'`), its counts against the trace formula at the 31 good primes below 150, and the conductor of `C'` |
| `stepH-thuemahler.m` | the branch fibre, exactly: the two branch points stay distinct modulo every prime of good reduction, which makes the weighted y-coordinate of the fibre a unit outside {2, 7, 11, 13} and the quintic factor of the star's sextic an S-unit — a Thue–Mahler equation whose 21 solutions (`stepH-thuemahler.sols`, from Gherga–Siksek's solver, 237 s) give 7 possible fibres, of which only the fibre at infinity passes the count test |

With these, *if* `Y` is hyperelliptic then it is `Y_2`, with nothing numerical in the chain. What
they do not do is prove hyperellipticity: that is the position of one rational divisor on the star,
and the note's Section 3 says which Borcherds form would settle it (the search is in progress).

The Thue–Mahler solver is not in this repository: it is `Code/TMSolver/solveThueMahler.m` of
`github.com/adelagherga/ThueMahler`, run as `solveThueMahler([8,1,90,-239,128,16], 1, [2,7,11,13] :
coprime := false)`; the `.sols` file is its output, and the script re-verifies every entry (its
completeness is the solver's theorem).

Two of the inputs need packages that are not in this repository: `Genus2Conductor` and
`endomorphisms`, both bundled in `edgarcosta/CHIMP`.  On lovelace that bundle is at
`~/GitHub/CHIMP`; its `MagmaPolred` submodule had to be checked out, and line 41 of
`MagmaPolred/MagmaPolred.m` lacked the semicolon after `end intrinsic`, without which nothing in
`endomorphisms` loads.
