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

## 4. What remains, and it is ONE LOCAL IDENTITY

The Weil representation factors over the orthogonal decomposition of the
discriminant form, and `(rho_A tensor rho_B)_{00} = (rho_A)_{00} (rho_B)_{00}`,
so the constant terms factor over primes. The two lattices of a support differ
at EXACTLY ONE prime -- at `p | D'` both are ramified; at `s` one is the
maximal split order and the other Eichler of level `s`. Hence

    m_s ct(theta_s) - m1 ct(theta_1) = (1/C) * prod_{p != s} f_p * g_s

    f_p(g) := (p-1) * ct_p^ram(g)             [ramified factor at p]
    g_p(g) := (p+1)/2 * ct_p^Eich(g) - 1      [level-raising at p]
    C      := 48 * 2^(w-1)

Dividing by `prod_{all p} f_p`, D'-independence of the left side is
EQUIVALENT to

        g_p(g) = lambda(g) * f_p(g),     lambda independent of p.          (**)

That is the theorem. (**) is established numerically by section 3 (the global
constancy holds to 120 digits, and the reduction above is exact algebra). What
is NOT yet done is deriving (**) in closed form from the local Whittaker /
local density formulas -- `EisensteinLocalFactors.m`, Kudla-Yang, with the
conventions pinned in `whittaker-sign-convention`. That is a finite,
prime-by-prime computation and it is the honest remaining gap.

## What this changes

Asking "which `D'` does the theory select?" is MALFORMED. The mass formula
makes every three-prime `D'` produce the same Eisenstein series, so
"odd part", "drop the smallest prime", "largest definite divisor", "larger
ramified prime" are not competing rules -- they are four names for one object,
and the telescoping is why.

The content is (*), and the mass-ratio law is its shadow in normalized
coordinates.
