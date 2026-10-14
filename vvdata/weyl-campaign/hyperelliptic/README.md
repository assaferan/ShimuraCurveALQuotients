# Hyperelliptic detection by Noether's count (2026-10-14)

A curve of genus g >= 3 is hyperelliptic iff the pairwise products of a basis of its holomorphic
differentials satisfy (g-1)(g-2)/2 linear relations (canonical image a rational normal curve);
otherwise there are exactly (g-2)(g-3)/2 (Max Noether). Two ways to get the differentials:

* a function field (`IsHyperellipticFunctionField`, HyperellipticDetection.m, PR #108): every curve
  stored in the CRV shape, and any cover whose double covers are stored (fibre product);
* for D = 1, the weight-2 cusp forms on Gamma_0(N) fixed by W, as q-expansions past the Sturm bound
  of weight 4 (`decide_d1.m`); controls `d1_controls.m`: Ogg's hyperelliptic X_0(40), (48), (47),
  (59) and the non-hyperelliptic X_0(34), (45), (64), (53), (55), all right.

Results on the 299 stored CRV curves (`hypsweep.m`, `hypclass.log`): 60 hyperelliptic, all by an
Atkin-Lehner involution; 58 over a conic without a rational point; the two covers of X_0(10,53) by
w_530 and w_53 over the rational conic y^2 = 5x^2 - 8x + 16 got models y^2 = f(x) (`hyp1053c.m`).

The 173 undecided curves of the 2026-09-30 rerun (`census_hyp.m`, `census_und*.m`): 7 at bases
with model files (all seven NOT hyperelliptic, `decide7.m`: 390_1 x2, 462_1, 510_1 x2, 22_15,
10_51); 7 at level 420 with the composite-level branch's files (all NOT hyperelliptic,
`decide420.log`: 6_35 x2, 10_21 x2, 14_15 x3); 22 at D = 1 (`decide_d1.log`: 19 NOT hyperelliptic,
3 hyperelliptic -- X_0(114)/{1,2,19,38}, X_0(130)/{1,2,13,26}, X_0(195)/{1,5,39,195}, the three
curves PR #51 found hyperelliptic with explicit models). The remaining 137: 132 at 114 bases with
D > 1 whose star curve has genus > 0 (no genus-0 quotient, so no Borcherds model and no fibre
product can exist), and 5 at genus-0 bases not yet run (15_14; 210_13 at M = 5460; 798_1, 870_1,
1326_1 at M = 1596, 1740, 2652).
