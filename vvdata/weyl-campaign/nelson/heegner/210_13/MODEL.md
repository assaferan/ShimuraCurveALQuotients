# X_0^210(13)/{1, 2, 21, 42, 65, 130, 1365, 2730}: y² = 4x⁸ − 20x⁷ + 97x⁶ − 20x⁵ − 138x⁴ + 20x³ + 97x² + 20x + 4

Built 2026-10-18 by the height route, TWICE, through two different elliptic quotients:

* via w₁₅: E = 2730s1, branch points of discriminants −312 and −120 over Q(√2), minus parts exactly
  2G′ (k = 2.000000 for both), one surviving configuration (b even), twist κ = −7555 unique with
  primes ≤ 800;
* via w₅: E = 2730u1 (2730u2 excluded by the lattice test, k = 2√2 there), branch points −420 and
  −1092 over Q(i), minus parts exactly 2G′, surviving configurations in the even class only, twist
  κ = 2 unique with primes ≤ 800.

**The two curves are Q-isomorphic** (`norm210.m`): different elliptic curves, different CM points,
different L-values, the same octic. The model above is the w₅ one after the coordinate change
x ↦ λx with λ = 2/15 and a shift (`small210.m`); it is invariant under x ↦ −1/x.

Checks that were predictions: bad primes exactly {2, 3, 5, 7, 13}; point counts agree with the trace
formula at the 104 primes 11 ≤ p ≤ 600; automorphism group of order 8; elliptic quotients 2730u1
(w₅), 2730s1 (w₁₅) and 210a2 (w₃ — a level-13 oldform, consistent with its fixed points of
discriminant −3640 lying over a quartic field, where the formula was not used). No involution has
rational fixed points, as the CM fields (Q(√2), Q(i), a quartic field) require.

Only b mod 2 matters: translating both branch points by a rational point of E gives an isomorphic
cover, and halvability depends on b mod 2; so the enumeration over b ∈ {0, 1} is complete and the
survivors at every even b in `pipe_210_w5.log` are translates of one curve. Rigorous given
Cai–Shu–Tian's theorem and the exact cover test. Hyperelliptic, as Nelson's numerics said.
