# PREREGISTRATION: the eleventh base X0^330(1)  (written BEFORE the measurement)

Date: 2026-08-26.  D = 330 = 2*3*5*11, N = 1, M = 660, |Dbar| = 217800 =
2^3*(3*5*11)^2, 2-part elementary (Z/2)^3.  First production run of the
Kubota closed-form rhov (eis32k.m, banked law b130990) -- the FFT could
never reach this base.

## What this base can and cannot decide

CANNOT separate "odd part" vs "largest definite divisor": for ANY even
squarefree D (even number of primes), the largest divisor with an odd
number of prime factors is D/2 = the odd part -- the two rules coincide
identically on even D.  They also agree with the unified reading
"q = D/p_min, s = p_min" which fits all ten banked bases.  The genuine
separator is an odd D with four primes (minimal D = 1155, M = 4620, level
4DN), currently blocked on (i) non-elementary 2-part (Z/4 Jordan
components -- extend Stromberg Def 2.15 handling in eis32k) and (ii) the
O(|D|)-per-YTab-key scans at |D| ~ 1e7 (need the CRT/Jordan closed form
for y0 and |D[c]|).

CANNOT separate scale laws tied to omega: 330 has the same profile as 210
(2 times three odd primes, omega(q) = 3), so "scale 2 universal" and
"scale 2^(omega(q)-2)" predict alike here.

WHAT IT IS: a REPLICATION test of the 210_1 re-split law at fresh primes
({3,5,11} vs {3,5,7}), the 11th consecutive base for the mass-ratio law,
and a second data point for the scale-2 reading.

## Preregistered predictions

P1 (support -- the re-split law): E_eis lies in span{ GROSS(165,1),
GROSS(165,2) } mod cusp; support {(165,1),(165,2)}, i.e. q = 165 = D/2,
s = 2, and NO third term (at N = 1 the (330,1) slot is indefinite --
formal mass 0, exactly as at 210_1).  Falsifiers: a fit needing
{(11,*)} (largest-prime rule), any {(33,*)}/{(55,*)} mixture, or a third
support term.  Gauge check as at f56c0a4: restricting the fit to the
competing support must FAIL rank.

P2 (weights): primary reading -- w = (-2, +3), i.e. 2x the s=2 prime law
(-1/(s-1), (s+1)/(2(s-1))) = (-1, 3/2), replicating the 210_1 scale 2.
Competing reading recorded in advance: constant |w|/mass = 8 as measured
at 210_1 (w1/mass1 = -8, w2/mass2 = +8 there), which here gives
w = (-10/3, +5).  The two SEPARATE at this base because the masses do
not scale by the same factor as the s=2-law weights.  Both readings
satisfy P3; the fit decides.

P3 (mass-ratio law, 11th base): |w2/w1| = mass(165,2)/mass(165,1) = 3/2.
Masses computed in advance (mass165.m, Eichler/Gross genus masses):
mass(165,1) = 5/12 (disc 27225, 4 classes), mass(165,2) = 5/8 (disc
108900, 4 classes).  Cross-check banked: (105,1) = 1/4, (105,2) = 3/8;
the 7 -> 11 prime swap multiplies both by (11-1)/(7-1) = 5/3 -- the
per-prime multiplicativity of the mass law.

P4 (structural, soft): kernelrat constants-map rank 15 and aliasing dim
25, as measured at 210_1 (the composite-D values; the prime-base
constants were rank 7 / aliasing 6).  Deviation here would be
informative, not fatal.

## Protocol

eis32k.m (closed-form rhov; the single a=0 word measured by FFT -- the
base-dependent S-coset sign bit is printed as AZERO) with
epool_330_1.txt (enum32m R=8 K=7, 169 quotients; character key = the
UNIVERSAL one verified on all seven banked mono dumps: s1 = s2 = 0 mod
24, parity bit 1 at p=2 only, via the synthetic carrier
mono330_1_synth.log).  Then kernelrat -> eisrat -> genallgross ->
allgross over all definite (D',R) with D'R | 330.  The verdict is the
rank-certified support + weights; SOLVE resid and "residual beyond
rank" are the arbiters (eisrat's RATIONAL banner alone is NOT trusted).

Validation gate run first: eis32k at 210_1 must reproduce the banked
FFT rhov (960/960 words incl. the S-coset via the AZERO-FFT path) and
the banked SOLVE resid ~ 1e-42 / identical BETAs.
