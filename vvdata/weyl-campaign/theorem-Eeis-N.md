# THEOREM (N > 1 prime): the three-term identity, proven — and lambda was a bug

Companion to `theorem-Eeis.md` (the `N = 1` theorem). Same method, same local
inputs; the new content is the tracked coset and the third theta term.

## Statement

Let `D = D' * s` be squarefree with exactly two prime factors (forced at
`N = 1` by indefiniteness of `B_D`, and the case of every banked `N > 1` base),
and let `N` be a prime with `gcd(D, N) = 1`. Let `L = L(D,N)` be the Shimura
curve lattice — ramified at `D'` and `s`, Eichler at `N` — and let `mu` be a
nonzero isotropic coset of `A_L`, which lives in the hyperbolic plane at `N`.
Write `Theta(D'',R) = mass(D'',R) * theta(D'',R)` with the **multiplicative**
mass

    mass(D'', R) = mass(D'', 1) * prod_{p | R, p prime} (p+1)/2 .

Then, at every coset,

    ct_L(gamma, mu)  =  1/2 * [Theta(D',Ns) - Theta(D',N)]
                             / (mass(D',Ns) - mass(D',N))
                        - 1/2 * theta(DN, 1) .

Nothing is fitted: both weights come out of the mass formula.

## First: lambda = 2N/(N+1) was NOT arithmetic

The banked `N > 1` residue — `lambda = 2N/(N+1)` holding on most cusp classes
and `1` on `{q,4,8,qr,4r,8r}` — was an artefact. The mass helper was called as
`mass([D'], [Rs])` with `Rs` the WHOLE Eichler level, contributing a single
factor `(Rs+1)/2` instead of the product over the PRIMES of `Rs`. At `N = 1`
the support is `(D',1)` and `(D',s)` with `s` PRIME, so the two agree and
twenty bases never saw it. At `N > 1` the level `Rs = N*s` is composite and
they diverge.

`lambda` could only ever have been a fudge: a single scalar cannot repair a
class-dependent error, which is exactly why it failed on a structured minority
of classes. `2N/(N+1)` and the class list are both artefacts and are
**withdrawn**. (`ctn1.m`, the first script of that sitting, passed `Rl`/`Rls`
as lists of primes and was right; `ctn1test.m` and `ctn1lambda.m` regressed.)

No `N = 1` result is affected — there `R = 1`, `Rs = s` prime.

## Proof

**(1) The factorisation is exact and the global factor is lattice-independent.**
`ctsplit.m` splits the closed form as

    ct_L(gamma) = glob(gamma) * prod_p loc_p^L(gamma) ,
    glob = u * (-a | |c|) * hilb(-a,c) ,
    loc_p = e(-sig8(p)/4) * [Gauss sum at p] * sqrt(dcn_p / n_p)   (xi2 at p=2)

`glob` does not see the lattice, so ratios and telescopings may be done prime by
prime. Self-check `glob * prod_p loc_p = ct` : **0 bad, worst 1.7e-119**, on all
four lattices at 22_3, 15_2, 10_7. (Valid for the zero coset and
`ord_2(c) != 1`, where the Stromberg `x_c` correction and the y-phase are
trivial; that is where the whole question lives.)

**(2) The tracked coset contributes 1 or 0 — no phase.** `mu` enters the closed
form only through the solvability of `c*y = mu` and the phase `e(a c Q(y))`;
`dcn` does not depend on it. The hyperbolic plane at `N` has exponent `N`, so
`mu in cA` iff `N` does not divide `c`; and then `y = c^{-1} mu` gives
`Q(y) = c^{-2} Q(mu) = 0` because `mu` is ISOTROPIC. Hence

    kappa := ct_L(gamma, mu) / ct_L(gamma, 0)  =  1  (N not| c),  0  (N | c).

Measured: **exactly `1.0000000`**, not merely modulus 1, at 22_3 / 15_2 / 10_7
(`ctn1local.log`).

**(3) The right-hand side telescopes to one product.** `G(D',N)` is Eichler at
`N` only; `G(D',Ns)` is Eichler at `N` and at `s`. With the multiplicative mass
their common factors `(N+1)/2` and `c_N^Eich` cancel out of the quotient, and
with `g_s = f_s/2` from the `N = 1` theorem (valid at every prime, `2`
included),

    T := [Theta(D',Ns) - Theta(D',N)] / (mass_Ns - mass_N)
       = c_{D'}^ram c_N^Eich * [ ((s+1)/2) c_s^Eich - 1 ] / ((s-1)/2)
       = c_{D'}^ram c_N^Eich * g_s / ((s-1)/2)
       = c_{D'}^ram c_N^Eich c_s^ram ,

while `theta(DN,1) = c_{D'}^ram c_s^ram c_N^ram`, all three primes being
ramified. So

    pred = 1/2 * c_{D'}^ram c_s^ram * ( c_N^Eich - c_N^ram ) .

**(4) That difference IS the tracked-coset indicator.** `c_N^Eich = 1` if
`N | c` else `1/N`; `c_N^ram = 1` if `N | c` else `-1/N`. Hence

    c_N^Eich - c_N^ram = 0 (N | c),   2/N (N not| c) ,

so `pred = 0` when `N | c`, and when `N` does not divide `c`

    pred = c_{D'}^ram c_s^ram / N = c_{D'}^ram c_s^ram c_N^Eich = ct_L(gamma,0) .

Comparing with (2): `ct_L(gamma,mu) = kappa * ct_L(gamma,0) = pred` in **both**
cases. QED.

The third term is not decoration: it is what makes the right-hand side vanish
on exactly the cosets the tracked coset kills.

**(5) The 2-parts.** As at `N = 1` the argument needs the four lattices'
2-adic normalisations to agree after telescoping. `ctsplit` shows it directly —
at 22_3, class 1: `L` and `G(66,1)` both give `loc_2 = (1+i)/4`, and the
telescoped 2-part of the `(11,·)` pair is `(1+i)/4` as well.

## Evidence

**Eight banked bases, every coset, no fitted scalar** (`ctn1exact.log`), and
the derived weights reproduce the banked empirical fits exactly:

    15_2 -1/2 T(5,2)   1   T(5,6)   -1/2 T(30,1)     288 cosets
    21_2 -1/2 T(7,2)   1   T(7,6)   -1/2 T(42,1)     384
    22_3 -1   T(11,3)  3/2 T(11,6)  -1/2 T(66,1)     576
    33_2 -1/2 T(11,2)  1   T(11,6)  -1/2 T(66,1)     576
    35_2 -1/4 T(7,2)   3/4 T(7,10)  -1/2 T(70,1)     576
    55_2 -1/4 T(11,2)  3/4 T(11,10) -1/2 T(110,1)    864
    77_2 -1/6 T(11,2)  2/3 T(11,14) -1/2 T(154,1)   1152
    10_7 -1   T(5,7)   3/2 T(5,14)  -1/2 T(70,1)     576

worst deviation 3.2e-119.

**PREDICTIVE TEST — eight bases never fitted, 16 support-instances, all
PREDICTED** (`ctn1predict.log`), worst 1.7e-119:

    6_5  10_3  14_3  6_7  14_5  26_3  39_2  34_3

with weights ranging over `-1/16 .. 3/2` (e.g. 34_3 at `D'=2`:
`-1/16 T(2,3) + 9/16 T(2,51) - 1/2 T(102,1)`) — all read off the mass, none
fitted. Control in the same run: dropping the third term fails at relative
1.41 to 3.50 at every one of the sixteen.

**Corollary, and it came true in that run: the support is NOT unique.**
`pred = c_{D'}^ram c_s^ram c_N^Eich` is symmetric in `D' <-> s`, so both
choices of `D'` give the same Eisenstein series. Both were tested at all eight
fresh bases and both hold — the `N > 1` analogue of the `210/330/1155`
finding ([[resplit-rule-is-not-a-selection]]).

**Corollary: the third term exists iff `omega(DN)` is odd.** With `omega(D) = 2`
and `N` prime, `omega(DN) = 3` always, so at `N > 1` in this family the third
term is always present — matching all twelve banked bases and confirming the
`N = 1` explanation (`B_D` indefinite forces `omega(D)` even).

## Scope, honestly

* `omega(D) = 2` and `N` PRIME. Not tested: `omega(D) = 4` with `N` prime, and
  `N` COMPOSITE.
* `N` composite is a genuinely different regime and the identity as stated
  **cannot** hold there: with two Eichler primes `N_1 N_2`, `omega(DN)` is even
  so `theta(DN,1)` does not exist, yet `kappa` still vanishes when `N_i | c`.
  The derivation says what to expect instead — a FOUR-term combination

      1/4 [ Tel(D'; N, Ns) - theta(D' N_1 s, N_2) - theta(D' N_2 s, N_1)
            + Tel(D' N_1 N_2; 1, s) ]

  since `prod_i (c_{N_i}^Eich - c_{N_i}^ram)/2` is `c_{N_1}^Eich c_{N_2}^Eich`
  off the bad cosets and `0` on them, and each of the four sign patterns is
  realised by a lattice with an ODD number of ramified primes. Untested.
* This is still the EISENSTEIN part only. `W(m) = a_{E*}(m)` modulo the cusp
  ideal remains the open piece ([[theorem-Eeis-odd-D]] item 3).
