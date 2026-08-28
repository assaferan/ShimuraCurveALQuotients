# Why there is no re-split rule: the derivation

Yesterday's finding was empirical -- every three-prime `D'` carries `E_eis`,
so the Gross-span test selects no support. This is *why*, and it turns the
observation into arithmetic.

## 1. The two banked laws are ONE identity

The campaign carried two separate empirical laws at `N = 1`: weight sum `1`,
and `|w2/w1| = rho` (the mass ratio). Writing `m1 = mass(D',1)`,
`m_s = mass(D',s)`, `rho = m_s/m1`, the fitted weights
`w1 = 1/(1-rho)`, `w2 = rho/(rho-1)` satisfy

    w1/m1 = -1/(m_s - m1)  and  w2/m_s = +1/(m_s - m1),   so  w1/m1 + w2/m_s = 0.

So in terms of the UNNORMALIZED mass-weighted genus sums `Theta = mass * theta`,

        E_eis  =  [ Theta(D',s) - Theta(D',1) ] / ( m_s - m1 ).            (*)

Two laws, one identity. The "weights" were the identity in normalized
coordinates.

## 2. The denominator is D'-independent -- PROVEN

    m1      = prod_{p|D'}(p-1) / (48 * 2^(w-1))
    m_s     = m1 * (s+1)/2
    m_s - m1 = m1*(s-1)/2 = prod_{p|D'}(p-1) * (s-1) / (2 * 48 * 2^(w-1))

`D'` and `s` PARTITION the primes of `DN`, so
`prod_{p|D'}(p-1) * (s-1) = prod_{p|DN}(p-1)`. It TELESCOPES: the value does
not depend on which prime is left out.

    DN = 1155: 480/384 = 5/4    at all four supports   [measured: 5/4]
    DN =  330:  80/384 = 5/24   at all four supports   [measured: 5/24]

For a ONE-prime `D'` there is no telescoping --
`m_s - m1 = m1*[prod_{q|R}((q+1)/2) - 1]` -- and indeed the values differ
(71/24, 47/12, 35/8, 115/24 at 1155). Those supports cannot join the identity,
which is exactly the one-prime/three-prime split observed.

## 3. The numerator is D'-independent -- proven at the S-cusp, verified everywhere

At the S-coset `ct(theta_L) = e(sig/8)/sqrt|A|` with `|A| = 2 (D'R)^2`, so

    m_s ct_s - m1 ct_1 = e(sig/8)/(D' sqrt2) * m1 * [ (s+1)/(2s) - 1 ]
                       = -e(sig/8) * prod_{p|D'}(p-1)(s-1) / (384 sqrt2 * D's)
                       = -e(sig/8) * prod_{p|DN}(p-1) / (384 sqrt2 * DN)

-- the same telescoping, and `D's = DN`. Predicted magnitudes 7.65267079206e-4
(1155) and 4.46405796203e-4 (330); both measured to 12 digits (`ctnumerator.m`).

Across ALL cosets the numerator agrees between supports to **1.8e-120**
(1155) and **2.6e-118** (330), relative.

## 4. The local identity, in closed form -- PROVEN

Constant terms factor over primes: the Weil rep is multiplicative on
orthogonal sums and `(rho_A ox rho_B)_00 = (rho_A)_00 (rho_B)_00`.  The two
lattices of a support differ at EXACTLY ONE prime (at `p | D'` both ramified;
at `s` maximal vs Eichler level `s`), so

    m_s ct_s - m1 ct_1 = (c_2/C) * prod_{p != s} f_p * g_s
    f_p := (p-1) c_p^ram ,   g_p := (p+1)/2 c_p^Eich - 1 ,   C := 48*2^(w-1)

**The local forms.**  At odd `p | DN`, both local discriminant groups are
`(Z/p)^2`, but the forms differ:

    Eichler level p (split, O = {c = 0 mod p}):  Q(v) = -x^2 - p yz
        => L*/L = (Z/p)^2 with  Qbar(b,c) = -bc/p      HYPERBOLIC
    Ramified (B = (u,p), u a nonsquare unit):    Nm(v) = -ux^2 - py^2 + upz^2
        => L*/L = (Z/p)^2 with  Qbar(b,c) = (uc^2 - b^2)/p   ANISOTROPIC

Verified by isotropic-vector counts on the actual Gross lattices: 1 for
ramified, `2p-1` for Eichler, at every prime and every support.

**The local Weil entries**, from the Gauss sums:

    hyperbolic:   sum_{b,c} e(-t bc/p) = p            => signature 0
    anisotropic:  sum e(t(uc^2-b^2)/p) = (u|p)(-1|p) g(1)^2 = -p
                  (u nonsquare so (u|p) = -1, and g(1)^2 = (-1|p) p)
                                                      => signature 4

so, with `c` the lower-left entry of the coset representative,

    c_p^Eich(g) =  1  if p | c ,   1/p  otherwise
    c_p^ram (g) =  1  if p | c ,  -1/p  otherwise

i.e. `c_p^ram = eps_p c_p^Eich` with `eps_p = +1` if `p | c`, `-1` otherwise.
Both verified against the lattices at 11520 cosets per prime, **0 mismatches**.

**Hence, in both cases at once,**

    g_p = (p+1)/2 c_p^Eich - 1 = (p-1)/2 * eps_p * c_p^Eich
    f_p = (p-1) c_p^ram        = (p-1)   * eps_p * c_p^Eich

        =>   g_p = (1/2) f_p ,   for EVERY p and EVERY coset.

(Measured: `g_p/f_p = [0.5, 0.5, 0.5, 0.5]` across `p = 3,5,7,11`.)

Therefore `prod_{p != s} f_p * g_s = (1/2) prod_{p | DN} f_p`, which does not
mention `s` at all.  **QED.**

## Why one-prime D' cannot join

For `D' = p` a single prime the partner is `(p, DN/p)`, Eichler at THREE
primes, so the bracket is `prod_q [(q+1)/2 c_q^Eich] - 1` -- a product of
three local factors minus one, which does NOT factor into `g_q`'s.  No
telescoping, no shared value, and indeed `m_s - m1` takes four different
values (71/24, 47/12, 35/8, 115/24 at 1155).  The one-prime/three-prime split
is thus derived, not observed.

## What this changes

Asking "which `D'` does the theory select?" is MALFORMED. The mass formula
makes every three-prime `D'` produce the same Eisenstein series, so
"odd part", "drop the smallest prime", "largest definite divisor", "larger
ramified prime" are not competing rules -- they are four names for one object,
and the telescoping is why.

The content is (*), and the mass-ratio law is its shadow in normalized
coordinates.
