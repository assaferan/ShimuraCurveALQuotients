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
