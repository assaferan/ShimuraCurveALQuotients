# THEOREM (odd D, N = 1): the identity is proven

Let `D` be odd squarefree with `w` prime factors (`w` even, since the Shimura
curve algebra `B_D` is indefinite), `N = 1`. Let `D'` be any product of `w-1`
of those primes and `s = D/D'`. Write `Theta(D',R) = mass(D',R) * theta(D',R)`
for the unnormalized mass-weighted genus sum. Then

        rho  =  [ Theta(D',s) - Theta(D',1) ] / ( mass(D',s) - mass(D',1) )

and in particular the right-hand side does not depend on the choice of `D'`.

This was an EMPIRICAL law fitted at twenty bases. It is now a theorem.

## Proof

**(1) Both sides factor over primes.** The Weil representation is
multiplicative on orthogonal sums and `(rho_A ox rho_B)_00 = (rho_A)_00
(rho_B)_00`, so for any lattice `M`, `ct(M)_g = prod_p c_p^M(g)` where
`c_p^M(g) := (rho_{A_{M,p}}(g)^{-1} e_0)_0`. And `rho` IS such an object for
the Shimura curve lattice at `N = 1`, where the tracked coset is the zero
coset.

**(2) The two local forms.** At odd `p`, both relevant local discriminant
groups are `(Z/p)^2`:

    Eichler level p (split):  O = {c = 0 mod p},  Q = -x^2 - p yz
        =>  Qbar(b,c) = -bc/p              HYPERBOLIC
    ramified:  B = (u,p), u a nonsquare unit,  Nm = -u x^2 - p y^2 + u p z^2
        =>  Qbar(b,c) = (u c^2 - b^2)/p    ANISOTROPIC

Verified on the actual lattices by isotropic counts (1 vs 2p-1).

**(3) The local Weil entries**, from the Gauss sums
`sum e(-tbc/p) = p` and `sum e(t(uc^2-b^2)/p) = (u|p)(-1|p) g(1)^2 = -p`
(using `(u|p) = -1` and `g(1)^2 = (-1|p)p`), i.e. signatures 0 and 4:

    c_p^Eich = 1 if p | c,  +1/p otherwise
    c_p^ram  = 1 if p | c,  -1/p otherwise            so  c_p^ram = eps_p c_p^Eich

Verified against the lattices at 11520 cosets per prime, 0 mismatches.

**(4) The local identity.** With `f_p := (p-1) c_p^ram` and
`g_p := (p+1)/2 c_p^Eich - 1`, in both cases at once

    g_p = (p-1)/2 eps_p c_p^Eich ,   f_p = (p-1) eps_p c_p^Eich   =>   g_p = f_p / 2 .

**(5) The right-hand side collapses.** The two lattices of a support differ at
exactly one prime, so with `C = 48*2^(w-2)`,

    mass_s ct_s - mass_1 ct_1 = (c_2^G / C) prod_{p != s} f_p * g_s
                              = (c_2^G / 2C) prod_{p|D} (p-1) prod_{p|D} c_p^ram
    mass_s - mass_1           = prod_{p|D} (p-1) / (2C)         [telescoping:
        prod_{p|D'}(p-1)*(s-1) = prod_{p|D}(p-1) since D' and s partition]
    =>  RHS = c_2^G * prod_{p|D} c_p^ram .

**(6) The left-hand side.** `B_D` is ramified exactly at the primes of `D`,
and ramification is a LOCAL condition -- it does not see the archimedean
place -- so the Shimura curve lattice has the SAME anisotropic plane at each
`p | D` as the definite Gross lattices do. Verified. Hence

        rho = c_2^L * prod_{p|D} c_p^ram .

**(7) The 2-parts agree.** For odd `D` both `A_{L,2}` and `A_{G,2}` are `Z/2`
with `4 Qbar = 1`, so `c_2^L = c_2^G`. Its Gauss sum is `1 + i`, normalized
`e(1/8)`, signature 1, giving `c_2 = 1` if `2 | c` and `e(1/8)/sqrt2`
otherwise -- exactly the two values observed.

Comparing (5) and (6): **rho = RHS**. QED.

## Numerical closure

At `D = 1155`: `c_2^L` and `c_2^G` computed independently by dividing the
known odd part out of `rho` and out of `ct(theta(385,1))` --
**13824 cosets, 0 mismatches, worst 1.45e-119**, and `c_2` takes exactly the
two predicted values `1` and `(1+i)/2`.

## What this does NOT cover

* **Even `D`.** Then `2 | D` is itself a ramified prime and step (3), derived
  for odd `p`, does not apply to it. The 2-parts still match at 330 (both
  `(Z/2)^3`, value multiset `[1,2,2,2,3,3,3]`), but a genuine 2-adic version
  of (3) is needed.
* **`N > 1`.** The tracked coset is no longer the zero coset and absolute
  normalisations shift (see the tracked-coset remark). The Eichler-at-`N`
  structure should enter the same local framework.
* **The Weyl vector itself.** This is the EISENSTEIN part. `W(m) = a_{E*}(m)`
  holds only modulo the cusp ideal, and aliasing is nonzero at some bases.
