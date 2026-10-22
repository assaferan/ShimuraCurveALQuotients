# Exact models of the genus-two Atkin–Lehner stars, from the periods

Produced by `periodstar.m` at 150 digits and certified by `certstar.m`. These are the stars under
the open genus-four quotients; the Nelson–Shimizu route could not reach them, because recognising
their invariants needs around a hundred digits and that route delivers twelve. The invariants below
have 20- to 31-digit numerators, which is the quantitative form of that gap.

## X_0^2002(1)\*

    y^2 + (x^2 + x) y = 16x^6 + 4x^5 + 180x^4 - 456x^3 + 196x^2 + 64x + 4

Cardona–Quer invariants, recognised with margins 2.9e92, 1.8e101 and 4.0e109 over the accidental
rate:

    4443541864293086393626096432993/8718043021712850681856
    -34439169873536221631669217/2177333422006206464
    186506522749409116433/543789565935616

Certified: #C(F_p) = p + 1 − a_p, untwisted, at all 13 primes in [3, 60] of good reduction. The
discriminant is supported on {2, 7, 11, 13}, exactly the primes of D = 2002 — which nothing in the
construction imposed.

## X_0^2030(1)\*

    y^2 + (x^2 + x) y = x^6 + 29x^5 + 13x^4 + 9x^3 + 159x^2 + 149x + 39

Invariants recognised with margins 1.4e95, 3.5e101 and 2.9e107:

    1200089628966417104818086601/3102846490240000000000
    123175783384857423892659/8865275686400000000
    63643057903538421401/25329359104000000

Certified at all 13 primes, untwisted, discriminant supported on {2, 5, 7, 29} = the primes of
D = 2030.

## Status of the other two

X_0^2418(1)\* and X_0^3390(1)\* are not done here. The runs reach the candidate and are then killed
in the certification step — see the warning below.

## ⚠ What the certification must not do

Computing Hecke operators up to p = 200 on the modular **abelian variety** at these levels is
expensive enough to get the process killed: the level-2002 run was terminated by SIGTERM
immediately after printing its candidate, losing the certification though the curve was correct.
Read the traces from the modular **symbols** factor instead, and stop at p = 60 — thirteen primes
already pin the twist and leave no freedom. `certstar.m` does this standalone, for a candidate
recovered from a log.

## Why these are only candidates until certified

The surfaces admit several principal polarizations — ten small combinations at each of these levels
— and different principal polarizations give different curves. A large recognition margin says the
invariants are genuinely rational, not that they belong to the star. Only the point counts against
the newform tie a candidate to the right object, and that is why no model here is quoted without
them.
