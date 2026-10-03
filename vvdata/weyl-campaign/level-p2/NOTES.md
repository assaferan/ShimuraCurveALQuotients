# The m = 0 analysis at a p^2-scaled lattice

What PR #66 proves at a level prime p || N, redone for p^2 || N -- the configuration the two
unobstructed non-squarefree bases (`6_25`, `6_49`) need, and the one Yang's Lemma 18 does not cover.

## The local lattice (p2local.m)

Built rather than assumed: the Eichler order of level N in B of discriminant D, its trace-zero
lattice, and the p-adic elementary divisors of the Gram matrix of tr(x ybar).

    p || N,  p odd   valuations 0, 1, 1      Yang's Lemma 18(1): unit line + p-scaled plane
    p^2 || N, p odd  valuations 0, 2, 2      unit line + p^2-SCALED plane
    4 || N           valuations 1, 2, 2      the usual extra power of 2 (tr vs Q)
    8 || N           valuations 1, 3, 3

⇒ **For p^e || N the negative plane is the p^e-scaled binary lattice.** That is the clean
generalisation of Lemma 18(1), and it is what the Whittaker computation needs.

## The local factor (Kudla-Yang Thm 4.3 with l_1 = l_2 = 2), CHECKED BY COUNTING (whit_p2.m)

At a p^2-scaled SPLIT plane (the case a firing discriminant forces, p split in k), the m = 0 factor
at the nonzero isotropic cosets takes THREE values, where X = p^(-s):

    cosets                                                 count        W_{0,p}        kappa^-_mu(0)
    exact order p^-2 (and isotropic)                       2p(p-1)      1              -log p / (p(p-1))
    exact order p^-1 with p | t_mu                         2(p-1)       1 + (p-1)X     -log p / (p-1)
    exact order p^-1 with p not| t_mu                      (p-1)^2      1 - X          0
    the zero coset                                         1            (1 + (p-2)X + (p-1)^2 X^2)/(1-X)

with t_mu = -Q(mu) summed over the non-integral components, and the total number of isotropic
cosets 3p^2 - 2p (65 at p = 5, confirmed by the count). The zero coset still has a SIMPLE POLE at
s = 0, as at level p, so the argument of prop:kappa0 goes through unchanged: with Phi(0) = -1,
kappa^-_mu(0) = -W_mu(1) log p / (p(p-1)), and W_0's value at X = 1 is p(p-1).

**Independently verified.** Every value above was also obtained by brute-force counting in the
normalisation of the preprint's Lemma lem:m0loc, W = (1-X)G(X) with
alpha_k = p^(-k) #{x in mu+L mod p^k L : Q(x) = 0 mod p^k}, at p = 3 and p = 5, one representative
per coset class (whit_p2.log). The zero coset's series 1, 4, 20, 20, 20, ... at p = 5 is exactly the
expansion of (1 + 3X + 16X^2)/(1-X), and the coset counts 24 = p^2-1, 8 = 2(p-1), 16 = (p-1)^2,
40 = 2p(p-1) all match.

## The consequence for the multiplier

Summing over the nonzero isotropic cosets, sum kappa^-_mu(0) = -4 log p, against -2 log p at level
p. The cosets with p not| t_mu drop out entirely. The orthogonal group of the discriminant form has
(at least) two orbits here -- exact order p^-1 and exact order p^-2 -- so the two surviving classes
can carry different constant terms, and with the -|CM(d)|/4 prefactor of Theorem B the correction
per CM point is

    mult = (1/2) ( c_1(0) + c_2(0) ),

c_1 the constant term on an order-p isotropic coset with p | t_mu, c_2 on an order-p^2 one. At
level p both collapse to the single class and this is the familiar (1/2) c_eta(0). ⚠ NOT yet
checked against any measured CM value: no base with p^2 || N has a model, which is the point of
pursuing `6_25` and `6_49`. The claim to test first is the coset-orbit count, i.e. that c_1 and c_2
are each constant on their class.

## The structure confirmed on the real lattices (isocount.m, isoclass.m, 2026-10-02)

The p^2-scaled-plane analysis above was derived from the model plane and checked by counting there.
It also makes a prediction about the ACTUAL lattices of the two bases that need it, testable in
seconds from `ShimuraCurveLattice` alone, with no Borcherds form in sight.

**Number of isotropic cosets of L^v/L.** 2N-1 at squarefree level; 3p^2 - 2p when p^2 || N, since
the D-part is anisotropic and the p-part is the p^2-scaled plane:

    6_5    |A| = 1800     isotropic 9     predicted 9     (2N-1)     MATCHES
    6_7    |A| = 3528     isotropic 13    predicted 13    (2N-1)     MATCHES
    15_2   |A| = 1800     isotropic 3     predicted 3     (2N-1)     MATCHES
    6_25   |A| = 45000    isotropic 65    predicted 65    (3p^2-2p)  MATCHES
    6_49   |A| = 172872   isotropic 133   predicted 133   (3p^2-2p)  MATCHES
    15_4   |A| = 7200     isotropic 10    predicted 8                MISMATCH

⚠ `15_4` differs because p = 2: its local valuations are (1,2,2), not (0,2,2), so the 2-adic count
is not 3p^2-2p. The formula is for ODD p -- which is the case 6_25 and 6_49 need. p = 2 at level 4
or 8 needs its own count.

**The finer split the multiplier rests on**, by the exact order of the coset at p:

    6_25 (p=5)   order p^0: 1    order p^1: 24    order p^2: 40     predicted 1, 24, 40
    6_49 (p=7)   order p^0: 1    order p^1: 48    order p^2: 84     predicted 1, 48, 84

matching 2(p-1) + (p-1)^2 = 24 resp. 48 of order p (of which 2(p-1) contribute and (p-1)^2 do not)
and 2p(p-1) = 40 resp. 84 of order p^2. So the two surviving classes are real and have the predicted
sizes on the actual lattices.

Also recorded: `ShimuraCurveLattice` works at both levels (probe625.m) with |A| = 2(DN)^2 as at
squarefree level, the star curves have genus 0 and the full curves genus 5 and 9, and the Hall
divisor group has 8 elements, {1,2,3,6,25,50,75,150} at 6_25.

## ⇒ WHY 6_25 IS OUT OF REACH: the pool, not the timeout (measured 2026-10-03)

The theory above is done, and nothing structural blocks the base: both lattices build, the
Atkin-Lehner group is already indexed by Hall divisors throughout the library, and the whole
quotient diagram of X_0^6(25) comes out instantly (16 curves, top genus 5, star genus 0). The
obstruction is the eta-quotient pool at level M = 300.

    pole order 55    1.6 s      299 lattice points      (a cached triple, reproduced exactly)
    pole order 135   ~45 min    4 075 460 lattice points   (272 MB of output)

and the form ring asks for fifteen orders, 135, 175, ..., 695. So the half-hour solver limit was
NOT the real obstacle -- with a four-hour budget the first deep solve finishes -- but it returns
FOUR MILLION generators, where the odd-D bases that do build have pools of seven or eight thousand.
The downstream echelon and kernel steps cannot take that, and the deeper orders are far worse.

⇒ **6_25 (and a fortiori 6_49, at M = 588) is out of reach by this route, for a structural reason
rather than a budget.** Do not raise NMZ_TIMEOUT and wait: the solve succeeding is what proves the
point. Level 300 has 18 divisors, so the polytope has 18 variables against 12 at M = 60, and the
pole orders needed are also an order of magnitude deeper; both push the point count up.

⇒ What this leaves. The level-p^2 m = 0 analysis above stands on its own and is verified, but it
cannot be tested against a measured CM value until some p^2-level base has a model, and none is
reachable with the present construction. The 15_4-style route is closed for the obstructed bases
(see ../nonsquarefree/NOTES.md) and the pool is the obstacle for the unobstructed ones. So the
non-squarefree frontier needs a different construction, not a bigger budget.

## N dividing the CONDUCTOR of d: a different plane, the same correction (2026-10-03)

Sachi's review of PR #66 points out that `SchoferFormula.m` applies the m = 0 term whenever N misses
the FUNDAMENTAL discriminant of d, while the proof assumes d itself fundamental. These differ when N
divides the conductor, and on X_0^15(2) that covers **d = -60, one of the three discriminants where
Guo-Yang's Table 45 validates the rule**. Reproduced and resolved:

**The lattice really is different** (`lminus.m`, counting isotropic cosets of L_- = L cap lambda^perp
straight from the lattice):

    d = -7, -15   fundamental        L_- at 2 is the EVEN plane 2*H     2 nonzero isotropic cosets
    d = -12, -60  conductor 2        L_- at 2 is 2*diag(u1,u2), u1 + u2 = 0 mod 4 (ODD type)
                                                                       1 nonzero isotropic coset

So the proof's count 2N-2 fails there, exactly as the review says.

**But the correction is the same** (`oddtype.m`, the same brute-force counting as above). At the odd
type 2*diag(1,-1):

    zero coset                 alpha_k = 1, 2, 3, 4, ...  so W_0 = 1/(1-X)   -- still a SIMPLE POLE
    the one isotropic coset    alpha_k = 2, 2, 2, ...     so W_nu = 1 + X    (value 2 at X = 1)
    the two anisotropic ones   alpha_k = 0                so W = 0

Then prop:kappa0's argument runs unchanged: W_nu/W_0 = 1 - X^2, which vanishes at s = 0 with
derivative 2 log 2, so **kappa^-_nu(0) = -2 log 2, exactly TWICE the -log N/(N-1) = -log 2 of the
fundamental case**, and with the -|CM(d)|/4 prefactor one coset at double weight gives

    -(1/4) * 1 * (-2 log 2) * c_nu(0)  =  (1/2) c_nu(0) log 2,

identical to the fundamental-d answer. ⇒ **The code is right at d = -60, and the proof extends to
N | conductor** rather than needing to exclude it; the paper should state this case.

⚠ **What this does NOT explain is the review's d = -12 disagreement** (for F = fs[-2] - fs[-1] the
code's log 2 part is off there). Since the local factor and the prefactor now match at -12 too, the
discrepancy must be in WHICH coset's constant term the code uses: `M0MultiplierExact` reads c(0) at a
nonzero isotropic coset of L^v/L, which does not depend on d, while the coset that actually occurs is
nu in L_-^v/L_-, and the correspondence between them DOES depend on d. At -60 the two evidently
agree and at -12 they need not. That is the next thing to check, and it is a question about the
multiplier's coset bookkeeping, not about the local Whittaker factor.

⚠ Note on the counting script: for an ANISOTROPIC coset alpha_0 = 0, not 1 (Q(mu) is not integral),
so the "1 +" in the generating function must be dropped there; with it one gets 1 - X instead of 0.
The isotropic rows, which are what matter, are unaffected.

## The inert conductor prime: NO pole at the zero coset, so NO m = 0 term (inertplane.m, 2026-10-03)

Counting on the ACTUAL negative plane `L_- = L ∩ λ^⊥` of `X_0^15(2)` (not on a guessed shape), the
zero coset's series `α_k` is

    d = -15   (fundamental, 2 split)   2, 3, 4, 5, 6, 7, 8     simple pole
    d = -60   (conductor 2, 2 split)   1, 2, 3, 4, 5, 6, 7     simple pole
    d = -240  (conductor 4, 2 split)   1, 2, 2, 4, 6, 8, 10    simple pole
    d = -12   (conductor 2, 2 INERT)   1, 2, 1, 2, 1, 2, 1     bounded: no pole
    d = -48   (conductor 4, 2 INERT)   1, 2, 2, 4, 2, 4, 2     bounded: no pole

The `m = 0` correction at a level prime comes from that pole (prop:kappa0), so at a level prime
dividing the conductor it fires iff the prime SPLITS in the CM field. Rule now in `SchoferFormula.m`.
Checked against Guo-Yang's Table 45 (fs[-2] value = (1280/9)|s(s-2)|): `-28, -60, -240` (fires) and
`-48` (does not) all exact. ⚠ The earlier conductor-4 defect (odd part halved) was separate: the
`m > 0` sum used `h(d)` of the order where the formula takes `h(d_0)` of the field (fixed, PR #66).

## Which cosets of the plane actually enter the formula at a conductor prime (fibre.m, 2026-10-03)

At a conductor prime `L_N` is NOT `L_+ (+) L_-`: measured index 1 (fundamental), 2 (f = 2), 8 (f = 4).
A coset `mu` of `L_-^v/L_-` enters the decomposition of some `phi_eta` only if `mu` lies in `L^v`.
Integral-norm cosets of the plane: `-15`: 2, both in `L^v` over nonzero `eta`; `-60`: 1, in `L^v`;
`-240`: 3, ONE in `L^v` (the others enter no `eta`); `-12`: 1, in `L^v`; `-48`: 3, ONE in `L^v`.
⇒ the pure `m = 0` piece at the inert plane is NOT zero in the formula (one coset over a nonzero
`eta`, with a nonzero first-order coefficient), so the observed absence of a `log 2` correction at
`-48` must come from the `x != 0` terms over the enlarged fibre (pairs with `mu_+ != 0`, which
exist because `L_+ = Z lambda` with `Q(lambda) = f^2 |d_0|`), not from the `m = 0` factor alone.
That accounting is the missing step of a proof.

## The lattice at a conductor prime, measured at p = 2 and p = 3 (conductor.m, 2026-10-03)

    base      d        d0    f   p at p   [L : L_+ (+) L_-]   L_- elem.div. p-vals   zero coset         integral cosets / in L^v
    15_2    -60      -15    2   split    2^1                 (1,1)                  1,2,3,4,5  POLE    1 / 1
    15_2    -240     -15    4   split    2^3                 (1,3)                  1,2,2,4,6  POLE    3 / 1
    15_2    -960     -15    8   split    2^5                 (1,5)                  1,2,2,4,4  POLE    7 / 1
    15_2    -12      -3     2   inert    2^1                 (1,1)                  1,2,1,2,1  none    1 / 1
    15_2    -48      -3     4   inert    2^3                 (1,3)                  1,2,2,4,2  none    3 / 1
    15_2    -192     -3     8   inert    2^5                 (1,5)                  1,2,2,4,4  none    7 / 1
    15_2    -28,-112 -7     2,4 split    2^1, 2^3            (1,1), (1,3)           POLE               1/1, 3/1
    10_3    -72      -8     3   split    3^1                 (0,2)                  1,3,5,7,9  POLE    2 / 2
    10_3    -648     -8     9   split    3^3                 (0,4)                  1,3,3,9,15 POLE    8 / 2
    10_3    -387     -43    3   inert    3^1                 (0,2)                  1,3,1,3,1  none    2 / 2
    10_3    -603     -67    3   inert    3^1                 (0,2)                  1,3,1,3,1  none    2 / 2

Pattern, now PROVED for odd p (standalone, lem:conductor): lambda = (p a', b; p^2 c', -p a') with b a
unit, L_-,p = <-1> ⊥ <-p^{2k}|d_0| v^2>, index p^{2k-1}, p^k - 1 integral cosets of which the p - 1
multiples of e_2'/p lie in L^v. Inert fundamental d at p || N: no optimal embedding (nu_p = 0), as
predicted; the inert case arises only through the conductor. p = 2: index 2^{2k-1}, divisors
(1, 2k-1), ONE integral coset in L^v -- computed, not yet derived.

## The restored terms at a conductor prime are a FIBRE SUM, and it is exact for all nine forms (fibresum.m, 2026-10-03)

The sum Theorem B restores, before the direct-sum identification of prop:mult, is

    T = sum_{nu != 0 in L_-^v/L_-}  kappa^-_nu(0)  sum_{x in L_+^v, x + nu in L^v}  c_{[x+nu]}(-Q(x)),      correction = -T/4 per point,

with kappa^-_nu(0) = log p * (W_nu/W_0)'(X = 1) from the counting series at p (prop:fibre in the
standalone).  `fibresum.m` evaluates it on X_0^15(2) for the nine forms of the model set; the local
series are reconstructed from alpha_0..alpha_10 (`fibresum.log`) and unchanged with alpha_0..alpha_12
(`fibresum12` run, scratch).  At a conductor prime two new kinds of pair appear: nu in L^v with
x = lambda_0/2 NOT in L but x + nu in L (coefficient c_0(-Q(lambda_0/2)) = c_oo(-|d|/16), the pole of
f at tau_{d/4}), and nu NOT in L^v with x = lambda_0/4, x + nu in L^v over a nonzero coset (the cusp-0
coefficient at q^{-3/4}).  With m = (1/2)c_eta(0), a = c_oo(-|d|/16), b = cusp-0 coefficient of q^{-3/4}:

    d = -240 (split, f = 4):  one coset in L^v, kappa = -2 log 2          -> (m + a) log 2
    d = -48  (inert, f = 4):  that coset kappa = -4/3 log 2, two cosets outside L^v kappa = -1/3 log 2
                                                                           -> (2/3 (m + a) + 1/3 b) log 2
    d = -60, -28 (split, f = 2): one coset, kappa = -2 log 2              -> m log 2  (= the fundamental recipe)

TRUTH for every form: each Borcherds form is a function on the genus-0 star curve with known divisor
(DivisorOfBorcherdsForm), so |value| = C prod |s(d) - s(d_i)|^{m_i} with s from Table 45 (27 points,
tests/_offline/GuoYang_15_2.m).  `fibrepipe.m` + the switches in `fibre-switches.diff` (applied to a
scratch copy of SchoferFormula.m, never to a worktree) give the raw table in variant A (current code)
and B (stripped at -240/-48 plus the fibre sum, injected DOUBLED -- `evidence.py` halves it); C is read
at d = -7 and checked at 6-8 further points per form before -240/-48 are judged (`evidence.log`):

    fibre sum:              18 / 18 exact (9 forms x {-240, -48})
    "fire iff split" rule:   8 / 18  (right at -240 iff a = 0, i.e. tau_-60 not in the divisor;
                                      right at -48 iff m + a = 0 and b = 0 -- INCLUDING the one form
                                      tests/M0PoleSum.m checks, which is why the rule looked right)
    Yang's conductor term (Kappa, `Yang_tt`): 2^(4/3) at -48 for the two forms with b != 0 (the
                                      pipeline's RationalNumber conversion then dies -- variant A crashed)

⚠ tests/M0PoleSum.m called its form `fs[-2]`; the table rows follow Keys(fs), and row 1 is key 11
(the cover W = {1,3,5,15}, divisor (-40)+(-120)-2(-12), hence (1280/9)|s(s-2)|).  Fixed on kappa0-proof.
⚠ The code still applies the splitting rule (exact at every conductor-2 point; conductor-4 points are
never offered to the model search).  Implementing T needs the local factors at a conductor prime in
closed form or by counting (brute force is fine at p = 2, 3; 9^k at p = 3 for k <= 8 or so).

## Closed form of the m = 0 local factors at an odd conductor prime (wcond_check.m, 2026-10-03)

On L_-,p = <-1> + <-p^{2k} c> (lem:conductor), with eps = (d_0/p) and nu_r = (0, r/p^k), rho = ord_p r:
alpha_j(nu_r) = p^{floor(j/2)} (j <= 2 rho), (1+eps) p^rho (j > 2 rho);
alpha_j(0) = p^{floor(j/2)} (j <= 2k), (1+eps)(p-1)p^{k-1} ceil((j-2k)/2) + p^{k-(j mod 2)} (j > 2k).
kappa^-_{nu_r}(0) = -2 p^{rho-k+1}/(p-1) log p (split), -2/(p^{k-1}(p+1)) (p^rho + 2(p^rho-1)/(p-1)) log p (inert).
`wcond_check.m` counts on the actual plane of X_0^10(3) at p = 3: d = -72, -648 (split, k = 1, 2),
-387, -603 (inert, k = 1), every 3-primary integral coset, 6 terms: all agree (`wcond_check.log`).
At p = 2 the same formulas give every series of `fibresum.log` (k = 1, 2, both types) -- lem:Wcond in
the standalone (proved for odd p; p = 2 computed).

## The 2-adic lattice at a conductor prime is DERIVED (standalone lem:conductor2, 2026-10-03)

lambda = (2a', 2b'; 4c', -2a') with b' a unit, c' even (after conjugating by w_2 -- the same step the
odd-p proof had skipped: optimality only says one of b, c is a unit), a'^2 + 2b'c' = -4^{k-1}|d_0|.
Plane: <-1> + <-4^{k-1} c>, c = |d_0|/b'^2 = 3 mod 4, -c = 1 mod 8 iff 2 splits.  Index
[L : Z lambda_0 (+) L_-] = 2^{2k-2} (DIRECT at conductor 2; the earlier "2^{2k-1}" used lambda, not
the primitive lambda_0).  Dual Z/2 x Z/2^{2k-1}; integral cosets (0, t/2^{k-1}) [rho = 1 + ord t] and
(1/2, t/2^k), t odd [rho = 0]: 2^k in all, the one of order 2 lies in L^v.  A direct Hensel count
(1, 2, 2(1+eps) roots of a unit square mod 2^n for n = 1, 2, >= 3) gives the odd-p series verbatim, so
lem:Wcond holds for every p.  This is what fibresum.log / inertplane.log counted.

## A conductor prime OUTSIDE the level: the term exists, and it is a pole sum over d/q^2, d/q^4, ... (qprime.m, 2026-10-03)

At odd q not dividing DN with q^k || f, L_q = M_2(Z_q) trace-zero is unimodular; after GL_2(Z_q)
conjugation lambda = (0, d; 1, 0), L_- = <-1> + <-|d|>, index q^{2k}, integral cosets (0, r/q^k) --
NONE in L^v -- and the fibre over nu_r is x = t lambda with t = -r/q^k mod Z_q: the CM vectors of
d/q^{2(k-rho)}, all over eta = 0 (standalone lem:unimod).  q not dividing d_0: lem:Wcond verbatim;
q | d_0: anisotropic, kappa^- = -2(q^{rho+1}-1)/((q-1)q^k) log q.  Implemented (M0FibreCorrection with
Unimodular := true, odd q only).

TEST PROBLEM: the nine-form divisor relations (qprime.m) are BLIND here -- the term is proportional to
the pole order at tau_{d/q^2}, i.e. linear in the divisor, so they hold with and without it.
What is NOT blind: the class polynomial.  At d = -588 (three star points, cubic field) the norms
N|s|, N|s-2|, N|(s+1/12)(s-5/4)| must come from a monic cubic H in Q[X] (classpoly.py): WITH the term
H = X^3 - 191/54 X^2 + 343/432 X - 83^2/(2^8 3^3 7), whose field has discriminant -588 and the SAME
reduced polynomial as FieldsOfDefinitionOfCMPoint gives (fod.m); WITHOUT it, no rational cubic.
Same at -1960 (field disc -1960).  7-adic Newton polygon: every root has ord_7 = -1/3 -- all three
points reduce to the pole tau_{-12}, the Kronecker congruence s(tau_{dq^2}) = s(tau_d) mod q from
the q-isogeny; at -1960 the roots are divisible by 7 (s(tau_{-40}) = 0), at -735 N(s - 5/4) acquires
the 7 (s(tau_{-15}) = 5/4).  ⚠ q = 2 outside an odd level (35_1, 39_1, 51_1, 55_1, 57_1, 87_1 tables
have such points: -12, -28, -48, -60) is NOT derived -- L_2 is not unimodular there.
GY tables with ODD such points: 21_2 (-100, q = 5), 55_1/58_1/82_1/94_1 (-27, q = 3), 87_1 (-147, q = 7).

## External check of the ramified branch: X_0^58(1), d = -27 (gyforce_58_1.m, 2026-10-03)

The GY tables never REACH their conductor-3/5/7 points (CandidateDiscriminants offers conductors 1
and 2 only), so `gyforce_58_1.m` forces d = -27 = 3^2 (-3) into the table of X_0^58(1) -- a base with
NO level prime, so the q = 3 term is the only m = 0 term -- and runs the cross-ratio check of
GuoYangCheck.m.  The term is log 3 on the Hauptmodul form (key -1) alone, and s(-27) matches
Guo-Yang's 3/19 (all 7 comparable discs OK in both rows).  3 | d_0: this is the anisotropic
(ramified) branch kappa^- = -2(q^{rho+1}-1)/((q-1)q^k) log q.
