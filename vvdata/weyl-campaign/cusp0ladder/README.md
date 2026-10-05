# The t-shift ladder at a fixed cusp-0 order (2026-10-05)

`ladder.m M NLO NHI MM W0` checks that the enumerated rung (M, NHI, MM) is spanned, as a space of
forms (rank of the q-expansions at oo), by the enumerated rung (M, NLO, MM) shifted by the weight-0
quotients with poles at oo only (`tshift_w0_M.txt`, from `nmzsolve.py M k 0 out 0 1 0`).  Results:
204 (20,32)->(45,32) 63/63; 60 (0|3|8,8)->(11,8) 17/17; 372 (32,60)->(92,60) 124/124; 380 (28,65)->(93,65) 132/132.
`fbrank.m` rank-tests nmzsolve.py's own fallback output against the truth (204: 274 points, 63/63;
with the m = 0 thinning it was 88 points, 61/63).
`span204*.m` are the FAILED routes: products of the m = 0 rung with weight-0 two-cusp quotients span
24/38 of (204,20,32) (span204c, span204e iterates); the Atkin-Lehner reversal of an m = 0 rung lands in
the character twisted at the odd primes (`nmzsolve_tw.py`, NMZ_ODDPAR=1 enumerates it:
`twisted_204_32_0.txt`), and even the right-character reversal plus the cut rung spans 14/38.
