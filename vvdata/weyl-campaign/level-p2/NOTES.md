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
