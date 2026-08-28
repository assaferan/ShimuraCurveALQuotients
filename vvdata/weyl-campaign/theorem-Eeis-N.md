# THEOREM (N > 1): the identity, proven — and lambda was a bug

Companion to `theorem-Eeis.md` (the `N = 1` theorem). Same method, same local
inputs; the new content is the tracked coset and the extra theta terms. The
general form below **subsumes the `N = 1` theorem** as its `S = {}` case.

## Statement (general form)

Let `D` be squarefree with `omega(D)` even, written `D = D' * s` with `s` any
one of its primes; let `N` be squarefree with `gcd(D,N) = 1`. Let `L = L(D,N)`
be the Shimura curve lattice — ramified at the primes of `D`, Eichler at the
primes of `N` — let `mu` be an isotropic coset of `A_L`, and put

    S = { p | N : mu_p != 0 } .

Write `Theta(Delta,R) = mass(Delta,R) * theta(Delta,R)` with the
**multiplicative** mass, and let `Tel_s` be the telescoped pair over `s`:

    mass(Delta, R) = mass(Delta, 1) * prod_{p | R, p prime} (p+1)/2

    Tel_s(Delta; R) = [Theta(Delta,Rs) - Theta(Delta,R)]
                      / (mass(Delta,Rs) - mass(Delta,R))

Then, at every coset of the level `4DN`,

    ct_L(gamma, mu)  =  2^(-|S|) * sum_{T subset S} (-1)^|T| Term(T) ,

    Term(T)  =  theta( D' s prod T , N / prod T )      |T| ODD
             =  Tel_s( D' prod T ; N / prod T )        |T| EVEN

Nothing is fitted: every weight comes out of the mass formula. `|T|` even is
exactly when `omega(D' s prod T)` is even — the algebra is INDEFINITE, there is
no Gross lattice — and that is precisely when telescoping over `s` is what
supplies the missing ramified factor.

Special cases:

* `N = 1`: the only isotropic coset is `0`, `S` is empty, and the formula reads
  `ct_L = Tel_s(D'; 1)` — the `N = 1` theorem of `theorem-Eeis.md`.
* `N` prime, `mu != 0`: `S = {N}`, and the two subsets give

      ct_L(gamma, mu) = 1/2 Tel_s(D'; N) - 1/2 theta(DN, 1) ,

  the banked three-term identity. That `theta(DN,1)` shape is special to `N`
  prime: in general the odd-`|T|` term is `theta(D' s prod T, N/prod T)`.

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
`dcn` does not depend on it. The hyperbolic plane at `p` has exponent `p`, so
`mu in cA` iff no `p in S` divides `c`; and then `y = c^{-1} mu` gives
`Q(y) = c^{-2} Q(mu) = 0` because `mu` is ISOTROPIC. Hence

    kappa := ct_L(gamma, mu) / ct_L(gamma, 0)  =  prod_{p in S} [ p not| c ] .

Measured: **exactly `1.0000000`**, not merely modulus 1, at 22_3 / 15_2 / 10_7
(`ctn1local.log`).

**(3) The Eichler-minus-ramified difference IS the tracked-coset indicator.**
By (1), `ct_L(gamma, 0) = prod_{p|D} c_p^ram * prod_{p|N} c_p^Eich`. The two
local values are `c_p^Eich = 1` if `p | c` else `1/p`, and `c_p^ram = 1` if
`p | c` else `-1/p`, so at every prime

    ( c_p^Eich - c_p^ram ) / 2  =  c_p^Eich * [ p not| c ] .

Multiplying that over `p in S` and using (2),

    ct_L(gamma, mu) = kappa * ct_L(gamma, 0)
      = prod_{p|D} c_p^ram * prod_{p|N, p notin S} c_p^Eich
        * prod_{p in S} ( c_p^Eich - c_p^ram ) / 2 .

**(4) Expanding that product over subsets gives exactly the right-hand side.**
The `T`-term of the expansion is
`prod_{p|D} c_p^ram * prod_{p in T} c_p^ram * prod_{p|N, p notin T} c_p^Eich`,
with sign `(-1)^|T|` and weight `2^-|S|`. It remains to see that this is
`Term(T)`:

* `|T|` ODD: `omega(D' s prod T) = omega(D) + |T|` is odd, so
  `G(D' s prod T, N/prod T)` exists — ramified at the primes of `D` and of `T`,
  Eichler at the rest of `N` — and its `ct` is exactly that product.
* `|T|` EVEN: that discriminant is indefinite. Instead take the pair
  `(D' prod T, N/prod T)` and `(D' prod T, (N/prod T) s)`, which differ at `s`
  alone. With the multiplicative mass all common factors cancel out of the
  quotient, and with `g_s = f_s/2` from the `N = 1` theorem (valid at every
  prime, `2` included),

      Tel_s(D' prod T; N/prod T)
        = [prod_{p|D'} c_p^ram prod_T c_p^ram prod_rest c_p^Eich]
          * [ ((s+1)/2) c_s^Eich - 1 ] / ((s-1)/2)
        = [ ... ] * g_s / ((s-1)/2)  =  [ ... ] * c_s^ram ,

  which restores the missing `c_s^ram` and gives the same product. QED.

So the extra theta terms are not decoration: they are what makes the right-hand
side vanish on exactly the cosets the tracked coset kills.

**(5) The 2-parts.** As at `N = 1` the argument needs the lattices' 2-adic
normalisations to agree after telescoping. `ctsplit` shows it directly — at
22_3, class 1: `L` and `G(66,1)` both give `loc_2 = (1+i)/4`, and the
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

**COMPOSITE `N` AND EVERY ISOTROPIC COSET — 112 of 112** (`ctn1general.log`).
The general form was then tested at `35_6` (`N = 2*3`) and `15_14`
(`N = 2*7`), both `s` choices, and — crucially — at **every** nonzero isotropic
coset, not only the convention's `iso[2]`: 14 and 38 of them respectively, plus
`22_3` as an `N`-prime regression. Worst 2.7e-119. The term count follows `|S|`
exactly as predicted:

    S = {2,3} at 35_6, s=7:  1/4 Tel(5;6,42) - 1/4 T(70,3)
                             - 1/4 T(105,2) + 1/4 Tel(30;1,7)
    S = {3}   at 35_6, s=7:  1/2 Tel(5;6,42) - 1/2 T(105,2)
    S = {2}   at 35_6, s=7:  1/2 Tel(5;6,42) - 1/2 T(70,3)

Note `#iso = prod_{p|N} (2p-1)` (15 at `N = 6`, 39 at `N = 14`), generalising
the banked `2N-1`.

**`omega(D) = 4` at `N > 1` — 4 of 4** (`ctn1omega4.log`). The smallest instance,
`D = 210`, `N = 11`, level **9240**, 27648 cosets, `|A| = 2.9e7`: both `s = 2`
(`D' = 105`) and `s = 3` (`D' = 70`), two tracked cosets each, worst 2.2e-118.
So the theorem's `omega(D)` is genuinely arbitrary even, not `2` in disguise.

## Scope, honestly

* `mu` is assumed ISOTROPIC. That is the only place `Q(c^{-1}mu) = 0` is used;
  an anisotropic tracked coset would pick up a genuine phase and is outside the
  statement.
* This is still the EISENSTEIN part only. `W(m) = a_{E*}(m)` modulo the cusp
  ideal remains the open piece ([[theorem-Eeis-odd-D]] item 3) — and it is now
  the ONLY remaining piece between the identity and the Weyl vector.
