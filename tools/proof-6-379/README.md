# The computation behind `paper/hyperelliptic-6-379.tex`

`X_0^6(379)*` is proved hyperelliptic by a Borcherds form on it, given by its principal part and
never constructed, whose divisor is a positive multiple of `P_-123 + P_-163 - P_-52 - P_-168`; the
Jacobian has no rational torsion, so the two degree-two divisors are linearly equivalent and the
curve has a `g^1_2`.

| file | what it does |
|---|---|
| `stepH-dual.m` | finds the principal part: computes the weight-3/2 cusp forms of level 4548 (377 of them, to precision 2402; about four hours, cached under `data/dual/` and not committed) and solves the linear system "pairs to zero with every cusp form, and the coefficients sum to the prescribed multiplicity at every discriminant carrying CM points" in the unknowns `c(-n)`, `n <= CAP`. Consistent at `CAP = 700`, not at `500`; the solution, cleared of denominators and made primitive, is `data/dual/cert_6_379.m` |

Run from the repository root:

    magma -b DD:=6 NN:=379 MM:=4548 BB:=2400 CAP:=700 TGT:="[<-123,1>,<-163,1>,<-52,-1>,<-168,-1>]" tools/proof-6-379/stepH-dual.m < /dev/null

The two tests check the certificate: `tests/DualCertificate6_379_divisor.m` (12 s) recomputes the
divisor of the lift from the principal part, checks that the four points are rational points of the
star, and recomputes the torsion bound from the trace formula; `tests/_offline/DualCertificate6_379.m`
(about an hour) recomputes the cusp forms and checks the pairing, and that a perturbed principal part
fails it.

Controls run with the same basis, not banked: the system is also consistent for the other three
fibre pairs the numerical screen found (`{-75,-228}`, `{-43,-312}`, `{-147,-408}`, each against
`{-52,-168}`), and it is not consistent for `P_-123 - P_-163`, two points of one fibre, whose class is
a nonzero point of a torsion-free group.
