// ============================================================================
// Table 10 (OanaFreddy): the 469 genus-2 curves X = X_0^D(N)/W that are NOT
// bielliptic via an Atkin-Lehner involution, and for which the authors were
// unsure whether X is GEOMETRICALLY bielliptic.
//
// Goal: for each (D,N,W) decide bielliptic-or-not heuristically.
//   1. Build the genus-2 model C = X_0^D(N)/W (same pipeline as ratpts_table6.m
//      / ratpts_table1.m: find the star curve for (D,N), EquationsOfCovers, then
//      pull out the cover whose group W = AllALsFromGens(gens, D*N)).
//   2. Take CHIMP's HeuristicDecomposition(C)[4][1] (the GEOMETRIC decomposition
//      descriptor: [dim,exp] pairs over Qbar) and flatten exp to a dimension multiset.
//        dims = [1,1]  => Jac splits geometrically into two elliptic factors
//                         => X is (geometrically) BIELLIPTIC.
//        dims = [2]    => Jac stays 2-dimensional (simple)
//                         => X is NOT bielliptic.
//   NB: the library intrinsic HeuristicDecompositionFactors uses the decomposition
//   over the BASE field Q (entry [3]), which misses splittings that only occur over
//   an extension -- so we read entry [4] (over the closure) directly instead.
//   (D=6,N=29 W=<w3,w29> is the user-anchored example; expected BIELLIPTIC.)
//
// METHOD CONSTRAINT: the equation pipeline requires N squarefree. The 40 rows of
// Table 10 with non-squarefree N are NOT in TABLE10 below (handled separately /
// out of scope here). 429 squarefree-N rows remain, grouped into 228 (D,N)
// star-curve groups; EquationsOfCovers is computed ONCE per group.
//
// RESOURCE ORDER (cheap -> expensive), per the project heuristics:
//   * larger N is harder; larger D*N is harder.
//   * monotonicity: fixing D and increasing N, or fixing N and increasing D, if
//     the first (D,N) is intractable the next almost surely is too -> stop early
//     along that ray and don't re-run.
//   TABLE10 is sorted by D*N ascending (then N, then D). Run from the top.
//
// USAGE (always wrap in an OS-level `timeout` on this shared machine):
//   magma idx:=15 ratpts_table10.m        // run ONLY group #15 (D=6,N=29)
//   magma lo:=1 hi:=20 ratpts_table10.m   // run groups #1..#20
//   magma maxdn:=400 ratpts_table10.m     // run every group with D*N <= 400
//   magma ratpts_table10.m                // no idx/range -> just prints the table
//
// LOGGING: append a one-line verdict per (D,N,W) to ratpts_table10_results.txt
// as runs land. Once a (D,N,W) has a recorded verdict OR a recorded
// reason-it-fails (timeout / not-enough-CM-points / LP-too-big), do NOT re-run
// it. See the RESULTS LOG block at the bottom of ratpts_table10_results.txt.
// ============================================================================

AttachSpec("ShimuraQuotients.spec");
AttachSpec("/Users/sachihashimoto/Repos/CHIMP/CHIMP.spec");
SetVerbose("ShimuraQuotients", 1);

LOGFILE := "ratpts_table10_results.txt";

// polymake LP dimension guard (ported from ratpts_table6.m -- SAME EquationsOfCovers
// pipeline, same feasibility wall). The Borcherds-form step enumerates lattice points
// of a polytope of dimension #Divisors(M), M = 4*D0 = 4*(D*N)/2^v2(D). Cost is driven
// by #Divisors(M), NOT by the pole-order bound n (BorcherdsForms.m's LP_SIZE_CUTOFF,
// default 10000, only bounds n). Empirically: #div<=12 completes reliably; #div=16-20
// only at small n; #div>=24 (3 odd prime factors of M) OOM-kills polymake (Killed:9)
// even at its minimum forced n -- and that kill takes down the whole magma process,
// aborting a sweep. So skip #div(M) >= DIV_CUTOFF up front: clean skip, not a crash.
// (For Table 10's 228 groups: 58 are #div<=12, 10 are #div 16-20, 160 are #div>=24.)
DIV_CUTOFF := 24;
polymake_level := func< D, N | 4 * ((D*N) div 2^Valuation(D, 2)) >;  // = M

// Each entry: <D, N, [ set-of-AL-subscripts generating W, ... ]>, D*N ascending.
TABLE10 := [*
    <57, 1, [ {3} ]>,  // DN=57
    <39, 2, [ {2,13} ]>,  // DN=78
    <26, 3, [ {3,13} ]>,  // DN=78
    <51, 2, [ {2,3} ]>,  // DN=102
    <106, 1, [ {2} ]>,  // DN=106
    <57, 2, [ {2,3}, {3,38}, {6,38} ]>,  // DN=114
    <38, 3, [ {2,3}, {3,19} ]>,  // DN=114
    <118, 1, [ {2} ]>,  // DN=118
    <26, 5, [ {2,5}, {10,13} ]>,  // DN=130
    <69, 2, [ {2,3} ]>,  // DN=138
    <77, 2, [ {2,11}, {2,77}, {14,22} ]>,  // DN=154
    <34, 5, [ {2,5}, {5,17} ]>,  // DN=170
    <10, 17, [ {10,34} ]>,  // DN=170
    <58, 3, [ {2,3} ]>,  // DN=174
    <6, 29, [ {3,29} ]>,  // DN=174  <-- user-anchored example, expect BIELLIPTIC
    <26, 7, [ {7,13} ]>,  // DN=182
    <14, 13, [ {2,13}, {7,13} ]>,  // DN=182
    <65, 3, [ {5,13}, {15,39} ]>,  // DN=195
    <39, 5, [ {3,13} ]>,  // DN=195
    <15, 13, [ {3,5}, {5,13} ]>,  // DN=195
    <35, 6, [ {2,5,7}, {6,7,10} ]>,  // DN=210
    <21, 10, [ {2,3,7}, {3,10,14} ]>,  // DN=210
    <15, 14, [ {2,3,5}, {5,6,14} ]>,  // DN=210
    <14, 15, [ {2,3,5}, {6,7,10} ]>,  // DN=210
    <6, 37, [ {6,74} ]>,  // DN=222
    <34, 7, [ {2,7}, {2,119} ]>,  // DN=238
    <82, 3, [ {2,41}, {3,82} ]>,  // DN=246
    <6, 41, [ {2,3}, {3,41} ]>,  // DN=246
    <86, 3, [ {2,43} ]>,  // DN=258
    <6, 43, [ {2,3}, {2,43} ]>,  // DN=258
    <91, 3, [ {3,7,13} ]>,  // DN=273
    <39, 7, [ {13,21} ]>,  // DN=273
    <21, 13, [ {3,13} ]>,  // DN=273
    <143, 2, [ {2,143}, {11,13}, {22,26} ]>,  // DN=286
    <22, 13, [ {2,11} ]>,  // DN=286
    <58, 5, [ {5,58} ]>,  // DN=290
    <302, 1, [ {2,151} ]>,  // DN=302
    <303, 1, [ {3,101} ]>,  // DN=303
    <14, 23, [ {7,23}, {7,46} ]>,  // DN=322
    <33, 10, [ {2,5,11}, {2,11,15}, {6,10,11} ]>,  // DN=330
    <22, 15, [ {2,3,11}, {2,11,15}, {3,5,11} ]>,  // DN=330
    <15, 22, [ {2,5,33}, {5,6,22} ]>,  // DN=330
    <10, 33, [ {2,3,5}, {2,3,11}, {3,10,22}, {6,10,22} ]>,  // DN=330
    <6, 55, [ {2,5,11}, {3,5,11} ]>,  // DN=330
    <339, 1, [ {3,113} ]>,  // DN=339
    <6, 59, [ {2,59} ]>,  // DN=354
    <6, 61, [ {2,3} ]>,  // DN=366
    <39, 10, [ {2,13,15}, {2,15,39}, {3,10,13}, {6,10,13} ]>,  // DN=390
    <26, 15, [ {2,3,13}, {2,5,39}, {3,5,13} ]>,  // DN=390
    <15, 26, [ {3,5,26}, {5,6,13}, {5,6,26} ]>,  // DN=390
    <10, 39, [ {2,3,5}, {2,3,13}, {5,6,26}, {6,10,13} ]>,  // DN=390
    <6, 65, [ {3,5,13} ]>,  // DN=390
    <57, 7, [ {3,133}, {7,57} ]>,  // DN=399
    <201, 2, [ {2,3,67} ]>,  // DN=402
    <203, 2, [ {2,7,29} ]>,  // DN=406
    <14, 29, [ {7,29}, {7,58} ]>,  // DN=406
    <62, 7, [ {2,217} ]>,  // DN=434
    <6, 73, [ {2,3}, {3,146} ]>,  // DN=438
    <26, 17, [ {2,221} ]>,  // DN=442
    <33, 14, [ {2,11,21}, {3,7,11} ]>,  // DN=462
    <22, 21, [ {2,3,7}, {2,3,77}, {2,7,11}, {3,7,11}, {3,7,22} ]>,  // DN=462
    <21, 22, [ {3,14,22} ]>,  // DN=462
    <6, 77, [ {2,3,7}, {2,7,11}, {2,11,21}, {3,7,11} ]>,  // DN=462
    <155, 3, [ {3,5,31} ]>,  // DN=465
    <235, 2, [ {2,5,47} ]>,  // DN=470
    <38, 13, [ {2,247} ]>,  // DN=494
    <51, 10, [ {2,3,85}, {3,5,17}, {3,10,34} ]>,  // DN=510
    <15, 34, [ {5,6,34} ]>,  // DN=510
    <10, 51, [ {2,3,17}, {5,6,34} ]>,  // DN=510
    <6, 85, [ {2,3,5}, {3,5,17} ]>,  // DN=510
    <514, 1, [ {2,257} ]>,  // DN=514
    <10, 53, [ {5,106} ]>,  // DN=530
    <546, 1, [ {2,3,7} ]>,  // DN=546
    <91, 6, [ {2,3,7,13} ]>,  // DN=546
    <26, 21, [ {3,13,14}, {6,13,14} ]>,  // DN=546
    <21, 26, [ {2,3,91}, {3,7,13}, {3,13,14} ]>,  // DN=546
    <14, 39, [ {2,7,13}, {2,13,21}, {3,14,26}, {6,7,13}, {6,7,26}, {6,14,26} ]>,  // DN=546
    <6, 91, [ {2,3,13}, {2,7,13}, {6,7,26} ]>,  // DN=546
    <33, 17, [ {3,11,17} ]>,  // DN=561
    <570, 1, [ {2,3,5}, {2,3,19}, {2,5,19} ]>,  // DN=570
    <10, 57, [ {2,5,19}, {3,10,19} ]>,  // DN=570
    <6, 95, [ {2,3,19}, {2,5,19}, {2,15,57}, {3,5,19}, {5,6,38}, {6,10,38} ]>,  // DN=570
    <14, 41, [ {2,7,41} ]>,  // DN=574
    <586, 1, [ {2,293} ]>,  // DN=586
    <614, 1, [ {2,307} ]>,  // DN=614
    <123, 5, [ {3,5,41} ]>,  // DN=615
    <214, 3, [ {2,3,107} ]>,  // DN=642
    <15, 43, [ {3,5,43} ]>,  // DN=645
    <327, 2, [ {2,3,109} ]>,  // DN=654
    <218, 3, [ {2,3,109} ]>,  // DN=654
    <14, 47, [ {2,7,47} ]>,  // DN=658
    <35, 19, [ {5,7,19} ]>,  // DN=665
    <134, 5, [ {2,5,67} ]>,  // DN=670
    <226, 3, [ {2,3,113} ]>,  // DN=678
    <690, 1, [ {2,3,5}, {5,6,46} ]>,  // DN=690
    <69, 10, [ {2,3,5,23} ]>,  // DN=690
    <10, 69, [ {2,5,23} ]>,  // DN=690
    <6, 115, [ {2,3,23}, {2,5,23}, {3,5,23}, {3,5,46}, {5,6,23} ]>,  // DN=690
    <142, 5, [ {2,5,71} ]>,  // DN=710
    <714, 1, [ {2,7,17} ]>,  // DN=714
    <119, 6, [ {2,3,7,17} ]>,  // DN=714
    <6, 119, [ {2,7,17}, {2,7,51}, {2,21,51}, {6,7,17}, {6,14,17} ]>,  // DN=714
    <65, 11, [ {5,11,13} ]>,  // DN=715
    <254, 3, [ {2,3,127} ]>,  // DN=762
    <770, 1, [ {7,10,22} ]>,  // DN=770
    <55, 14, [ {2,5,7,11} ]>,  // DN=770
    <35, 22, [ {2,5,7,11} ]>,  // DN=770
    <22, 35, [ {2,5,77}, {5,14,22}, {7,10,11} ]>,  // DN=770
    <10, 77, [ {2,11,35}, {10,11,14} ]>,  // DN=770
    <111, 7, [ {3,7,37} ]>,  // DN=777
    <798, 1, [ {3,7,19} ]>,  // DN=798
    <14, 57, [ {2,3,133} ]>,  // DN=798
    <26, 31, [ {2,13,31} ]>,  // DN=806
    <10, 83, [ {2,5,83} ]>,  // DN=830
    <858, 1, [ {2,11,13}, {3,11,13} ]>,  // DN=858
    <39, 22, [ {2,3,11,13} ]>,  // DN=858
    <26, 33, [ {2,3,11,13} ]>,  // DN=858
    <6, 143, [ {2,11,13}, {2,11,39}, {2,33,39}, {3,11,26}, {3,22,26} ]>,  // DN=858
    <870, 1, [ {2,3,5} ]>,  // DN=870
    <6, 145, [ {2,3,29}, {2,15,29}, {3,5,58}, {6,10,58} ]>,  // DN=870
    <46, 19, [ {2,19,23} ]>,  // DN=874
    <15, 59, [ {3,5,59} ]>,  // DN=885
    <910, 1, [ {2,7,13}, {2,13,35}, {7,10,13} ]>,  // DN=910
    <65, 14, [ {2,5,7,13} ]>,  // DN=910
    <26, 35, [ {2,5,7,13} ]>,  // DN=910
    <930, 1, [ {2,15,93}, {3,5,31}, {3,10,62}, {5,6,31} ]>,  // DN=930
    <62, 15, [ {2,3,5,31} ]>,  // DN=930
    <15, 62, [ {2,3,5,31} ]>,  // DN=930
    <10, 93, [ {2,3,155}, {3,5,62}, {5,6,62} ]>,  // DN=930
    <6, 155, [ {3,10,31}, {3,10,62}, {5,6,62} ]>,  // DN=930
    <6, 157, [ {2,3,157} ]>,  // DN=942
    <966, 1, [ {2,3,23}, {6,7,23} ]>,  // DN=966
    <6, 161, [ {2,21,23}, {6,7,23}, {6,14,23}, {6,14,46} ]>,  // DN=966
    <326, 3, [ {2,3,163} ]>,  // DN=978
    <58, 17, [ {2,17,29} ]>,  // DN=986
    <14, 71, [ {2,7,71} ]>,  // DN=994
    <10, 101, [ {2,5,101} ]>,  // DN=1010
    <146, 7, [ {2,7,73} ]>,  // DN=1022
    <22, 47, [ {2,11,47} ]>,  // DN=1034
    <6, 173, [ {2,3,173} ]>,  // DN=1038
    <1110, 1, [ {2,3,37}, {2,5,37}, {3,10,74} ]>,  // DN=1110
    <1122, 1, [ {2,3,17}, {2,11,17}, {6,11,34} ]>,  // DN=1122
    <33, 35, [ {3,5,7,11} ]>,  // DN=1155
    <1190, 1, [ {2,17,35} ]>,  // DN=1190
    <34, 35, [ {2,5,7,17} ]>,  // DN=1190
    <1218, 1, [ {2,21,87} ]>,  // DN=1218
    <58, 21, [ {2,3,7,29} ]>,  // DN=1218
    <6, 203, [ {2,7,87}, {3,14,29}, {6,7,58} ]>,  // DN=1218
    <1230, 1, [ {2,3,41}, {2,5,41} ]>,  // DN=1230
    <82, 15, [ {2,3,5,41} ]>,  // DN=1230
    <1254, 1, [ {6,19,22} ]>,  // DN=1254
    <38, 33, [ {2,3,11,19} ]>,  // DN=1254
    <1290, 1, [ {2,15,129}, {3,10,86}, {5,6,86} ]>,  // DN=1290
    <10, 129, [ {2,3,5,43} ]>,  // DN=1290
    <1302, 1, [ {2,3,31}, {2,21,31}, {3,7,31}, {3,7,62}, {6,14,31} ]>,  // DN=1302
    <62, 21, [ {2,3,7,31} ]>,  // DN=1302
    <51, 26, [ {2,3,13,17} ]>,  // DN=1326
    <34, 39, [ {2,3,13,17} ]>,  // DN=1326
    <38, 35, [ {2,5,7,19} ]>,  // DN=1330
    <1410, 1, [ {6,10,94} ]>,  // DN=1410
    <6, 235, [ {2,5,141} ]>,  // DN=1410
    <26, 55, [ {2,5,11,13} ]>,  // DN=1430
    <6, 239, [ {2,3,239} ]>,  // DN=1434
    <26, 57, [ {2,3,13,19} ]>,  // DN=1482
    <6, 251, [ {2,3,251} ]>,  // DN=1506
    <22, 69, [ {2,3,11,23} ]>,  // DN=1518
    <14, 109, [ {2,7,109} ]>,  // DN=1526
    <14, 111, [ {2,3,7,37} ]>,  // DN=1554
    <6, 263, [ {2,3,263} ]>,  // DN=1578
    <1590, 1, [ {2,15,53}, {2,15,159}, {3,10,53}, {6,10,106} ]>,  // DN=1590
    <46, 35, [ {2,5,7,23} ]>,  // DN=1610
    <14, 115, [ {2,5,7,23} ]>,  // DN=1610
    <1722, 1, [ {2,21,123}, {6,7,82} ]>,  // DN=1722
    <1770, 1, [ {3,10,59} ]>,  // DN=1770
    <10, 177, [ {2,3,5,59} ]>,  // DN=1770
    <1794, 1, [ {2,23,39}, {6,13,23}, {6,13,46} ]>,  // DN=1794
    <14, 129, [ {2,3,7,43} ]>,  // DN=1806
    <1870, 1, [ {2,5,11,17} ]>,  // DN=1870
    <1914, 1, [ {2,3,319} ]>,  // DN=1914
    <6, 319, [ {2,3,11,29} ]>,  // DN=1914
    <1938, 1, [ {3,17,38}, {6,19,34} ]>,  // DN=1938
    <1974, 1, [ {2,3,329}, {2,21,47}, {3,7,94} ]>,  // DN=1974
    <6, 329, [ {2,3,7,47} ]>,  // DN=1974
    <2010, 1, [ {2,3,335}, {5,6,67} ]>,  // DN=2010
    <6, 335, [ {2,3,5,67} ]>,  // DN=2010
    <2130, 1, [ {2,5,213}, {2,15,71}, {3,10,71}, {3,10,142} ]>,  // DN=2130
    <10, 221, [ {2,5,13,17} ]>,  // DN=2210
    <770, 3, [ {2,7,11,15}, {2,11,15,21}, {3,5,7,22}, {5,6,7,11}, {5,6,7,22} ]>,  // DN=2310
    <462, 5, [ {3,5,7,22}, {3,5,14,22}, {3,7,10,11}, {5,6,7,11}, {5,6,11,14}, {5,6,14,22} ]>,  // DN=2310
    <330, 7, [ {2,3,5,7}, {2,3,7,11}, {2,3,35,55}, {2,5,21,33}, {2,7,15,33}, {3,5,7,22}, {3,10,14,22}, {6,7,10,11}, {6,10,11,14} ]>,  // DN=2310
    <210, 11, [ {2,3,5,11}, {2,3,35,55}, {3,5,7,11}, {3,5,11,14}, {3,7,10,11}, {6,10,11,14}, {6,10,14,22} ]>,  // DN=2310
    <14, 165, [ {2,3,5,7,11} ]>,  // DN=2310
    <6, 391, [ {2,3,17,23} ]>,  // DN=2346
    <1365, 2, [ {2,3,5,7,13} ]>,  // DN=2730
    <546, 5, [ {2,3,7,65}, {2,3,13,35}, {2,3,35,65}, {2,13,15,21}, {2,15,21,39}, {5,6,7,13} ]>,  // DN=2730
    <390, 7, [ {2,3,5,91}, {2,3,35,65}, {2,5,13,21}, {3,10,14,26}, {5,6,7,26}, {6,7,10,13}, {6,7,10,26}, {6,10,13,14}, {6,10,14,26} ]>,  // DN=2730
    <210, 13, [ {2,5,7,39}, {2,5,21,39}, {2,13,15,21}, {3,10,13,14}, {5,6,7,13}, {6,7,10,26}, {6,10,13,14} ]>,  // DN=2730
    <10, 273, [ {2,3,5,7,13} ]>,  // DN=2730
    <6, 473, [ {2,3,11,43} ]>,  // DN=2838
    <6, 497, [ {2,3,7,71} ]>,  // DN=2982
    <3210, 1, [ {2,3,5,107} ]>,  // DN=3210
    <3458, 1, [ {2,7,13,19} ]>,  // DN=3458
    <3486, 1, [ {2,3,7,83} ]>,  // DN=3486
    <3542, 1, [ {2,7,11,23} ]>,  // DN=3542
    <510, 7, [ {2,3,17,35}, {2,3,35,85}, {2,5,17,21}, {3,5,7,34}, {3,5,14,34}, {5,6,14,17}, {6,10,14,17} ]>,  // DN=3570
    <210, 17, [ {2,3,7,85}, {2,5,7,51}, {2,7,15,17}, {3,5,14,17}, {5,6,7,34}, {6,10,14,17} ]>,  // DN=3570
    <3930, 1, [ {2,3,5,131} ]>,  // DN=3930
    <1330, 3, [ {2,3,5,7,19} ]>,  // DN=3990
    <570, 7, [ {2,3,5,133}, {2,5,7,57}, {2,15,19,21}, {3,5,7,38}, {3,7,10,19}, {3,10,14,38}, {5,6,7,19}, {6,7,10,19}, {6,10,14,38} ]>,  // DN=3990
    <210, 19, [ {2,3,7,95}, {2,5,7,57}, {2,5,19,21}, {2,7,15,57}, {3,10,14,38}, {6,10,14,38} ]>,  // DN=3990
    <6, 665, [ {2,3,5,7,19} ]>,  // DN=3990
    <4002, 1, [ {2,3,23,29} ]>,  // DN=4002
    <4278, 1, [ {2,3,23,31} ]>,  // DN=4278
    <330, 13, [ {2,3,5,143}, {2,11,13,15}, {3,5,22,26}, {5,6,11,13} ]>,  // DN=4290
    <4326, 1, [ {2,3,7,103} ]>,  // DN=4326
    <4710, 1, [ {2,3,5,157} ]>,  // DN=4710
    <4890, 1, [ {2,3,5,163} ]>,  // DN=4890
    <1122, 5, [ {2,3,5,11,17} ]>,  // DN=5610
    <858, 7, [ {2,3,7,11,13} ]>,  // DN=6006
    <870, 7, [ {2,3,5,7,29} ]>,  // DN=6090
    <1254, 5, [ {2,3,5,11,19} ]>,  // DN=6270
    <930, 7, [ {2,3,5,7,31} ]>,  // DN=6510
    <390, 17, [ {2,3,5,13,17} ]>,  // DN=6630
    <690, 11, [ {2,3,5,11,23} ]>,  // DN=7590
    <1110, 7, [ {2,3,5,7,37} ]>,  // DN=7770
    <714, 11, [ {2,3,7,11,17} ]>,  // DN=7854
    <462, 17, [ {2,3,7,11,17} ]>,  // DN=7854
    <210, 43, [ {2,3,5,7,43} ]>   // DN=9030
*];

procedure log_row(D, N, gens, verdict, model)
    fh := Open(LOGFILE, "a");
    fprintf fh, "%o\t%o\t%o\t%o\ty^2=%o\n", D, N, gens, verdict, model;
    delete fh;  // flush/close
end procedure;

// ---------------------------------------------------------------------------
// Classify one genus-2 model via CHIMP's heuristic Jacobian decomposition.
// ---------------------------------------------------------------------------
procedure check_group(C, gens, D, N)
    desc := Sprintf("D=%o N=%o W=<%o>", D, N, gens);
    g := Genus(C);
    if g ne 2 then
        printf "  [%o] WARNING genus = %o (expected 2); skipping bielliptic test\n", desc, g;
        log_row(D, N, gens, Sprintf("SKIP-genus=%o", g), "n/a");
        return;
    end if;
    f := HyperellipticPolynomials(C);
    printf "  [%o] model y^2 = %o\n", desc, f;
    // CHIMP HeuristicDecomposition(C) returns [* Kiso, Kdecinfo, base-dec, geo-dec *].
    // Entry [4][1] is the GEOMETRIC descriptor: a list of [dim, exp] pairs, one per
    // distinct simple factor of Jac over Qbar (dim = #Rows = factor dimension, exp =
    // multiplicity). Geometric biellipticity needs the decomposition over the CLOSURE
    // (entry [4]), not over the base field Q (which is what the library
    // HeuristicDecompositionFactors uses via entry [3]). Flatten exp to get the
    // dimension multiset: [[1,1],[1,1]] and [[1,2]] both -> [1,1] (E1xE2 / E^2);
    // [[2,1]] -> [2] (geometrically simple).
    decgeodesc := HeuristicDecomposition(C)[4][1];
    dims := Sort(&cat[ [ pair[1] : i in [1..pair[2]] ] : pair in decgeodesc ]);
    if dims eq [1,1] then
        verdict := "BIELLIPTIC";
        printf "  [%o] Jac splits (dims %o) => %o\n", desc, dims, verdict;
    elif dims eq [2] then
        verdict := "NOT-bielliptic";
        printf "  [%o] Jac simple (dims %o) => %o\n", desc, dims, verdict;
    else
        verdict := Sprintf("UNCLEAR(dims=%o)", dims);
        printf "  [%o] unexpected decomposition dims %o\n", desc, dims;
    end if;
    log_row(D, N, gens, verdict, f);
end procedure;

// ---------------------------------------------------------------------------
// Run one (D,N) group: compute cover equations once, then test each W.
// ---------------------------------------------------------------------------
procedure run_group(entry, curves)
    D := entry[1]; N := entry[2]; gensets := entry[3];
    printf "\n==== D=%o N=%o (D*N=%o) ====\n", D, N, D*N;
    if not IsSquarefree(N) then
        printf "  N=%o is not squarefree; method N/A; skipping\n", N;
        return;
    end if;
    M := polymake_level(D, N);
    ndiv := #Divisors(M);
    if ndiv ge DIV_CUTOFF then
        printf "  polymake level M=%o has #div=%o >= %o; OOM-doomed, skipping\n",
            M, ndiv, DIV_CUTOFF;
        return;
    end if;
    t0 := Realtime();
    if not exists(Xstar){X : X in curves | X`D eq D and X`N eq N and IsStarCurve(X)} then
        printf "  no star curve found for (D,N)=(%o,%o); skipping\n", D, N;
        for gens in gensets do log_row(D, N, gens, "SKIP-no-star-curve", "n/a"); end for;
        return;
    end if;
    try
        tgts := { AllALsFromGens(gens, D*N) : gens in gensets };
        // Cheap predict-and-skip: if X* lacks the CM points the targets need (e.g. a
        // genus-2 target needs 2g+5=9 but only 4 exist), no amount of Borcherds work
        // can determine the equation. Skip in ~tens of seconds instead of ~15 min.
        enough, need, have := EnoughCMPointsForTargets(Xstar, curves, tgts);
        if not enough then
            printf "  insufficient CM points (need=%o, have=%o); skipping before Borcherds work\n", need, have;
            reason := Sprintf("SKIP-insufficient-CM(need=%o,have=%o)", need, have);
            for gens in gensets do log_row(D, N, gens, reason, "n/a"); end for;
        else
            crv_list, ws, keys := EquationsOfCovers(Xstar, curves : Targets := tgts);
            printf "  computed %o cover equations (targets-restricted) in %o s\n", #crv_list, Realtime()-t0;
            for gens in gensets do
                W := AllALsFromGens(gens, D*N);
                if not exists(k){k : k in keys | curves[k]`W eq W} then
                    printf "  [D=%o N=%o W=<%o>] not among computed covers (keys); skipping\n", D, N, gens;
                    log_row(D, N, gens, "SKIP-not-in-covers", "n/a");
                    continue;
                end if;
                idx := Index(keys, k);
                C := crv_list[idx];
                check_group(C, gens, D, N);  // logs its own verdict/skip row
            end for;
        end if;
    catch e
        printf "  ERROR on (D,N)=(%o,%o): %o\n", D, N, e`Object;
        reason := Sprintf("FAILED:%o", e`Object);
        for gens in gensets do log_row(D, N, gens, reason, "n/a"); end for;
    end try;
    printf "  ---- group (D=%o,N=%o) done in %o s ----\n", D, N, Realtime()-t0;
end procedure;

// ---------------------------------------------------------------------------
// Entry point: pick what to run via command-line assignments.
// ---------------------------------------------------------------------------
printf "Table 10 driver: %o (D,N) groups, sorted by D*N ascending.\n", #TABLE10;

if assigned idx then
    i := StringToInteger(idx);
    curves := GetHyperellipticCandidates();
    printf "Loaded %o candidate curves. Running single group #%o.\n", #curves, i;
    run_group(TABLE10[i], curves);
elif assigned lo then
    a := StringToInteger(lo);
    b := assigned hi select StringToInteger(hi) else #TABLE10;
    curves := GetHyperellipticCandidates();
    printf "Loaded %o candidate curves. Running groups #%o..#%o.\n", #curves, a, b;
    for i in [a..b] do
        run_group(TABLE10[i], curves);
    end for;
elif assigned maxdn then
    cap := StringToInteger(maxdn);
    curves := GetHyperellipticCandidates();
    printf "Loaded %o candidate curves. Running groups with D*N <= %o.\n", #curves, cap;
    for entry in TABLE10 do
        if entry[1]*entry[2] le cap then
            run_group(entry, curves);
        end if;
    end for;
elif assigned rest then
    curves := GetHyperellipticCandidates();
    // Complement of the default reprioritized run: everything we did NOT do by default,
    // i.e. N in {2,3} or N composite. Sorted D*N ascending so the fast #div/insufficient-CM
    // skips come first. With the DIV_CUTOFF + point-cap + CM pre-check guards this is safe.
    comp := [ e : e in TABLE10 | not (#PrimeDivisors(e[2]) eq 1 and e[2] notin {2, 3}) ];
    Sort(~comp, func< a, b | (a[1]*a[2]) ne (b[1]*b[2]) select (a[1]*a[2]) - (b[1]*b[2])
                             else (a[2] ne b[2] select a[2] - b[2] else a[1] - b[1]) >);
    printf "Loaded %o candidate curves. Running the REST: %o group(s) not in the default reprioritized set (N in {2,3} or composite N); D*N ascending.\n",
        #curves, #comp;
    for entry in comp do
        run_group(entry, curves);
    end for;
else
    curves := GetHyperellipticCandidates();
    // Reprioritized run order: keep only N with a SINGLE prime factor (N prime, since N
    // is squarefree here) and drop N in {2,3}. Then sort small-N-first (D ascending as a
    // stable tiebreaker within an N).
    prio := [ e : e in TABLE10 | #PrimeDivisors(e[2]) eq 1 and e[2] notin {2, 3} ];
    Sort(~prio, func< a, b | a[2] ne b[2] select a[2] - b[2] else a[1] - b[1] >);
    printf "Loaded %o candidate curves. Reprioritized to %o group(s) (single-prime N, N notin {2,3}); small N first.\n",
        #curves, #prio;
    for entry in prio do
        run_group(entry, curves);
    end for;
end if;

printf "\nDONE.\n";
exit;
