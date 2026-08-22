# The theta_g phase law from Stromberg's explicit formula — derivation notes

Target (measured, cuspclass2.py, 570 checks / 6 bases): at the canonical class-g coset
(bottom row (g,1)) of Gamma_0(M), for odd g | D coprime to N,

    (rho(w^{-1}) e_0)_{eta*} = (g/sqrt|A|) e(theta_g),
    8 theta_g = 1 + 6 r(g) + 4 s(g) mod 8,
    r(g) = #{p|g : p = 3 mod 4},  s(g) = #{p<q|g : (p/q)(q/p) = +1};

dead iff N|g; at 4|g real with sign mu(g/4), magnitude g/M.

Backbone: Stromberg, Math. Z. 275 (2013) 509–527, Thm 6.4 + Def 6.1 + Lemmas 3.6–3.9,
6.10, 6.11; metaplectic cocycle via Thm 4.1 (Hilbert symbols). Our machinery's rep
(VVWeilFFT Dual := true) is Stromberg's rho-tilde for Dbar = (A, -q): the code's
S-action constant e(-1/8)/sqrt|A| equals sigma_w(Dbar) = e_8(-sign(Dbar)) because
sign(Dbar) = -sign(A) = +1 mod 8 for these signature (1,2) lattices. [rho5.m verifies
the identification coset-by-coset, up to the metaplectic word unit.]

## Established steps

**1. Single-term collapse.** Thm 6.4 with x = 0:
    rho(A~) e_0 = xi(a,c) sqrt(|D[c]|/|D|) sum_{y in D/D[c]} e(ac Q(y)) e(B(x_c,y)) e_{cy + x_c}.
The eta*-component is a SINGLE term: nonzero iff eta* in cD + x_c, with y_0 any fiber
element.

**2. Support law.** For N | c: the N-part of cD is trivial (c kills the level-N
hyperbolic plane, where eta* has its nonzero component) and x_c is 2-adic (Def 2.15),
so eta* not in cD + x_c: entry = 0. This is the "class dead iff N | g" law.

**3. Magnitude law.** |entry| = sqrt(|D[c]|/|D|).
   * odd g: |D[c]| = prod_{p|g} p^2 (the odd Jordan pairs killed by c), so
     |entry| = g/sqrt|D| = g/sqrt|A|.  [15_2, g=5: sqrt(25/1800) = 5/(30 sqrt2). ok]
   * 4|g: the 2-part (Z/2)^3 is killed too: |D[c]| = 8 prod_{p | g odd} p^2, and
     sqrt(8 (g/4)^2 / (M^2/2)) = g/M.  [ok on 10_7 g=4,20; 6_11 g=4,12.]

**4. Isotropy collapse of the fiber phase.** For odd c, x_c = 0 and y_0 may be chosen
inside the N-part with c y_0 = eta*; then Q(y_0) = c^{-2} Q(eta*) = 0 (eta* isotropic),
so e(ac Q(y_0)) = 1 and the WHOLE phase is xi(a, c).

**5. xi(a,c) at the canonical rep** (Def 6.1; w = (a_w, b_w; g, 1), A = w^{-1} has
a = 1, c = -g). xi = e_4(-sign(Dbar)) xi_0 xi_2 prod_J xi(J):
   * e_4(-sign(Dbar)) = e_4(-1) = -i  (base-independent).
   * xi_2 = 1 (c odd).
   * xi_0 = ((-a)/c)(-a, c)_infty = ((-1)/(-g)) (-1,-g)_infty = (-1)^{(g-1)/2}
     = (-1)^{r(g)}  (odd squarefree g).
   * p | D odd, p NOT dividing c: the Jordan pair contributes
     G(c,0;A_p^{t1}) G(c,0;A_p^{t2}) = (t1 t2 / p) e_8(2-2p) = -1 for EVERY p
     (anisotropic: (t1 t2/p) = -(-1/p); e_8(2-2p) = (-1)^{(p-1)/2}); for p | N
     (hyperbolic pair) the same computation gives +1.
   * p | g: xi(J) = ((-a)/|J|) G(-ac, x_c; J) = 1 per component (the Gauss sum
     G(-ac, 0; A_p^t) with p | ac has q_c = p: trivial symbol, e_8(0)).
   * 2-part: three order-2 components A_2^{t_i}: G(c, 0; A_2^t) = (tc/2) e_8(tc)
     (Lemma 3.7, q = 2, c odd): depends on c = -g mod 8 and the t_i (fixed by the
     genus of the 2-part - same for all our even-DN bases).

So at the canonical rep,
    xi(1, -g) = [-i] * [(-1)^{r(g)}] * [(-1)^{#{p | D odd, p not dividing g}}]
                * [prod_i (t_i(-g)/2) e_8(t_i(-g))].
The D-dependent factor (-1)^{omega_odd(D) - omega(g)} must cancel against the
METAPLECTIC WORD UNIT u(w) (VVRhoInvE0FFT computes the word lift, not the canonical
lift; u = product of Kubota cocycles, Hilbert-symbol-valued by Thm 4.1) and/or
normalization — because the measured theta_g is base-independent. The s(g)
reciprocity term must also live in u(w): with a = 1 there is no Jacobi-symbol content
left in xi, and s(g) is NOT a function of g mod 8, so no product of the local factors
above can produce it. rho5.m's word-unit distribution (per class) is the direct
measurement of u(w).

## The 15_2 decode (rho5.m canonical-rep breakdown; all args in 1/8 turns)

    g   meas   u   eps   g3   gm   gd     [meas = u+eps+g3+gm+gd, all four classes OK]
    1    1     0    4    -3   -1    1
    3   -1     4    0    -3   -3    1
    5    1     4    4    -3    3    1
    15   3     4    0    -3    1    1

Identified closed forms (canonical rep (a,b,c,d) = (1, *, -g, 1); m = c d'^{-1} mod lev,
m = c mod 8; d' = cn - d = 1 mod 8):
* gd = G(d',0;Dbar) = (d'/|D|) e_8(sign(Dbar)) = e_8(1) ALWAYS: (d'/2) = +1 since
  d' = 1 mod 8, and (d'/(DN)^2) = 1.  [PROVEN]
* g3 = e_8(-3 sign(Dbar)) = e_8(-3).  [PROVEN]
* eps = ((-1)/(-g)) (-1,-g)_infty = -(-1/g): arg8 = 4 + 4r(g).  [PROVEN]
* gm = G(m,0;Dbar) factors over Jordan components: p | g gives 1; p | D odd, p not | g
  gives -1 (anisotropic pair, m-independent); p | N gives +1 (hyperbolic pair); the
  2-part gives G2(-g) with arg8 = -1 if -g = 3 mod 4 else +1 (genus-fixed function of
  -g mod 8; extraction from the table).  So arg8(gm) = 4(omega_D - omega(g)) + arg G2(-g),
  omega_D = #{p | D odd}.  [PROVEN modulo the G2 table]
* u (Kubota word unit of the canonical VVSTWord): 15_2 measured (0,4,4,4). For the
  total to be base-independent u MUST carry 4 omega_D mod 8. PREREGISTERED for 10_7
  (omega_D = 1): u(g=1) = 4, u(g=5) = 0 -- OPPOSITE of 15_2's g=1.
The s(g) term necessarily lives in u (with a = 1 nothing else can produce it; s is not
a function of g mod 8).

## Remaining steps
6. Pin u(w) for the canonical words: compute the cocycle product along VVSTWord of the
   canonical rep, reduce via Thm 4.1 / Lemma 4.3 to Hilbert symbols; show
   u(w) * (base-dependent factors of xi) = base-independent * e_8(4 s(g) + ...).
7. Assemble 8 theta_g = 1 + 6r + 4s and cross-check against the nine measured classes.
8. The 4|g case: x_c nonzero (2-adic), Lemma 6.10/3.7 casework; derive mu(g/4) sign.
9. Class-constancy (contrib constant): multiplier balance — e_0 is a
   rho(Gamma_0(M))-eigenvector by Lemma 5.5/5.6 (chi_D character) matching the
   eta-multiplier of f; T-invariance of the eta*-component from isotropy.

## RESOLUTION (same sitting): the law is DERIVED

10_7 breakdown (287/287 cosets pass, zero failures): u(g=1) = 0 there too -- the
omega_D-compensation happens in the 2-part genus, not in u:
* G2b (2|D bases: 10_7, 6_11) = -G2a (2|N bases: 15_2, 21_2, 33_2, 35_2), and
* omega_odd(D) is even iff 2 not| D (indefinite quaternion: omega(D) EVEN always),
so 4(omega_D) + argG2 is the SAME function on every base: gm-arg = -4 omega(g) +
argG2a(-g).  Base-independence PROVEN (not assumed).

**THE DERIVED LAW (canonical metaplectic lift, lower row (-g, 1)):**

    (rho_Dbar(A~) e_0)_{eta*} = (g / sqrt|A|) e_8( E0(g) ),
    E0(g) = 2 + 4 r(g) - 4 omega(g) + argG2a(-g)  mod 8,
    argG2a(x) = x + 4 [ (2/x) = -1 ]  mod 8      [G2a(x) = (2/x) e_8(x)].

The word-based machinery (VVRhoInvE0FFT on VVSTWord) differs by the metaplectic
center u = 4 [g > 1] (measured on all 6 canonical words of 15_2/10_7; cancels
against the a0-side branch in contrib, so the class-assembly value is convention-free).

**Verification: E0 + u reproduces ALL NINE measured class phases**
(g = 1,3,5,7,11,15,21,33,35 across six bases).

**A new falsifiable prediction:** at omega(g) >= 3 the derived law DIVERGES from the
fitted 1 + 6r + 4s (which was only ever constrained by omega <= 2 data): e.g. for g
with three prime factors all = 1 mod 4, derived 1 vs fitted 5. The preprint's
phase table should be restated via E0; "1 + 6r + 4s" is its omega <= 2 shadow
(equality checked for ALL class-multisets with omega <= 2, failure from omega = 3).

## Remaining to make it a theorem end-to-end
(i) G2a from the 2-adic Jordan type of the quaternion discriminant form in the two
    parities (Lemmas 3.7/3.8 casework) -- currently fitted on 7 residue points
    across 3 bases with the G2b = -G2a flip.
(ii) u = 4[g>1] for the canonical words via the Kubota cocycle (Thm 4.1) -- only
    needed for the word convention; the invariant statement needs nothing.
(iii) The 4|g classes (x_c nonzero, Lemma 6.10): same program, gives mu(g/4) g/M.
(iv) The 15_2 c = 2 mod 4 coset mismatches in rho5 (20 cosets; dead classes, no
    effect on the law) -- likely the x_c convention at 2||c; 10_7 has no failures.

## Step (i) CLOSED (g2jordan.m, ten bases): the 2-adic genera identified
The 2-part of Dbar is one of exactly TWO Jordan types, constant within parity:
* 2|N (15_2, 21_2, 33_2, 35_2):  Q-multiset {0,0,0, 1/4 x3, 1/2, 3/4} = A_2^1 + C_2,
  whence G2(c) = [(c/2)e_8(c)] * [1] = (2/c) e_8(c)          (Lemmas 3.7, 3.8);
* 2|D (10_7, 6_11, 6_5, 22_3, 34_5, 14_3): {0, 1/4, 1/2 x3, 3/4 x3} = A_2^1 + B_2,
  whence G2(c) = [(c/2)e_8(c)] * [(3/2)] = -(2/c) e_8(c).
Measured G2 phases match on every base (args (1,-1,1,-1) vs (-3,3,-3,3) at c=1,3,5,7).
The genus flip G2b = -G2a is therefore DERIVED, and with the quaternion parity theorem
the base-independence of E0 is complete at lemma level.

## Step (iii) CLOSED: the 4|g classes derive the same way
For 4 | c the level-4 component A_2^1 is killed (4 Q(y) integral on the 2-part), so
x_c = 0 again and the fiber phase collapses as before; magnitude sqrt(|D[c]|/|D|) = g/M
(established). The phase at the canonical rep (a = 1, c = -g):
  e_4(-1) [-2]  +  xi_0 = ((-1)/(-g))(-1,-g)_infty  +  xi_2 = -2(c_2 - 1 + s_2),
  c_2 = -g/4 (signed odd part), s_2 = sign of the 2-part (1 for A+C, 5 for A+B)
  +  4 #{p | D odd, p not | g}   [anisotropic pairs; hyperbolic give 0; p | g give 0;
                                  the 2-part Gauss factors are 1 by Lemmas 3.7/3.8 at 4|c].
All terms are EVEN: the entry is real, as measured. Verified against all four measured
4|g classes: 10_7 g=4 (arg 0 = mu(1)), g=20 (arg 4 = mu(5)); 6_11 g=4 (arg 0), g=12
(arg 4 = mu(3)); in each case with the word unit u = 4 (u = 4[g>1] holds for 4|g too).
The mu(g/4) sign is the anisotropic-pair count plus xi_0/xi_2 odd-part bookkeeping.

STATUS: the ENTIRE rho-side class law (support, magnitudes, odd-g phase E0, 4|g sign)
is now derived from Stromberg's formula, with only the Kubota word unit u = 4[g>1]
(convention-only; cancels in contrib) and the c = 2 mod 4 x_c convention (dead classes)
left as bookkeeping.

## The constancy character CONFIRMED (rho6.m)
For every gamma in Gamma_0(lev) tested (16 matrices incl. -I, T-conjugates, k = +-1, 2),
rho~(gamma) e_0 has off-diagonal mass < 1e-59: e_0 IS an eigenvector, with chi(gamma) in
{+-1, +-i} recorded per gamma. Class-constancy of contrib = [chi = eta-multiplier of f,
forced by well-definedness of the Guo-Yang average] + [T-invariance at eta* from isotropy].
The Kubota unit is definitional (Thm 4.1); its canonical-word value 4[g>1] is a finite
continued-fraction computation (still open as explicit casework, convention-only).
