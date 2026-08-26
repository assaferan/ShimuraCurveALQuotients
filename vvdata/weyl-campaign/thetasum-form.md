# The theta-sum closed form of the s-law (2026-08-26)

From massprobe.log (nine bases) + the Siegel-Weil confirmation
(genusfit15sw): over mass-weighted theta SUMS (not averages) the s-law is

    E_eis = c * [ -Theta(q,N) + Theta(q,Ns) - (N+1)/(N-1) * Theta(qsN,1) ]
            mod cusp,     c = 96 / ((s-1)(q-1)(N+1))

with the empirical mass identities (all nine bases, Gross-genus = Eichler/4):

    mass(q,N)    = (q-1)(N+1)/96
    mass(q,Ns)   = (q-1)(N+1)/96 * prod_{p|s}(p+1)/2
    mass(qsN,1)  = (q-1)(s-1)(N-1)/192

Consistency: w1 = -c*m1 = -1/(s-1); w2 = +c*m2 = (s+1)/(2(s-1)) at prime s
(generalizing to prod_{p|s}(p+1)/2^omega(s) / (s-1) at composite s);
w3 = -c*m3*(N+1)/(N-1) = -1/2 IDENTICALLY (the (N-1) cancels).

## Sharpened Law C preregistration for 210_1 (q=7, s=30, N=1)

c = 96/(29*6*2) = 8/29; mass(7,1) = 1/8; mass(7,30) = 9/8 (mass210.log).

    w1 = -1/29        (coincides with Law A -- both carry 1/(s-1))
    w2 = +9/29        (Law A says 31/58; Law B says 9/2)
    NO third Gross term: mass(qsN,1) has the factor (N-1) = 0, and indeed
    (210,1) is an indefinite structure (four primes) -- the genus does not
    exist. Law A's "-1/2 on (210,1)" is incoherent at N=1.

Discriminants: |w2/w1| = 9 (C) vs 31/2 (A) vs 9/8 (B).
