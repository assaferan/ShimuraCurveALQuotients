// The 34 Atkin-Lehner quotients X_0^D(N)/W of genus at least three that the filters of the
// pipeline leave undecided: no stage proves or rules out hyperellipticity for them.
//
// Provenance: the list of paper/hyperelliptic-candidates.tex (the theorem on the 34 and Table
// tab:candidates): the open quotients on which the differentials satisfy Noether's count for a
// hyperelliptic curve.  33 of them are hyperelliptic; X_0^141(2)/{1,6,94,141} is not.  The same 34
// are the 61 left open by the rerun of 2026-10-14 minus the 27 settled by the metric dual graph (22),
// a non-hyperelliptic component of a special fibre at p || N (3) and the models of X_0^6(35) (2).
//
// Entries are <D, N, W> with W the full Atkin-Lehner group, as a sorted sequence of the indices m.
// Read with eval (Read("data/open_34.m") cat "\nreturn open_34, open_gaps;").
open_34 := [
    <141, 2, [1, 3, 94, 282]>,
    <141, 2, [1, 6, 94, 141]>,
    <55, 7, [1, 35, 55, 77]>,
    <39, 11, [1, 11, 39, 429]>,
    <15, 38, [1, 2, 15, 19, 30, 38, 285, 570]>,
    <57, 10, [1, 5, 6, 19, 30, 95, 114, 570]>,
    <15, 46, [1, 5, 6, 23, 30, 115, 138, 690]>,
    <21, 34, [1, 3, 14, 17, 42, 51, 238, 714]>,
    <77, 10, [1, 2, 5, 7, 10, 11, 14, 22, 35, 55, 70, 77, 110, 154, 385, 770]>,
    <6, 185, [1, 2, 15, 30, 111, 185, 222, 370]>,
    <6, 187, [1, 3, 11, 17, 33, 51, 187, 561]>,
    <22, 51, [1, 6, 17, 22, 33, 102, 374, 561]>,
    <21, 58, [1, 2, 3, 6, 7, 14, 21, 29, 42, 58, 87, 174, 203, 406, 609, 1218]>,
    <39, 34, [1, 2, 3, 6, 13, 17, 26, 34, 39, 51, 78, 102, 221, 442, 663, 1326]>,
    <1326, 1, [1, 17, 78, 1326]>,
    <46, 31, [1, 2, 23, 31, 46, 62, 713, 1426]>,
    <6, 247, [1, 2, 19, 38, 39, 78, 741, 1482]>,
    <6, 253, [1, 2, 23, 33, 46, 66, 759, 1518]>,
    <6, 301, [1, 3, 14, 42, 43, 129, 602, 1806]>,
    <6, 323, [1, 6, 19, 34, 51, 114, 646, 969]>,
    <2002, 1, [1, 2, 13, 26, 77, 154, 1001, 2002]>,
    <2030, 1, [1, 10, 14, 35, 58, 145, 203, 2030]>,
    <6, 341, [1, 3, 11, 33, 62, 186, 682, 2046]>,
    <6, 379, [1, 2, 3, 6, 379, 758, 1137, 2274]>,
    <6, 385, [1, 5, 6, 11, 14, 21, 30, 55, 66, 70, 105, 154, 231, 330, 770, 1155]>,
    <6, 403, [1, 2, 13, 26, 93, 186, 1209, 2418]>,
    <2418, 1, [1, 2, 13, 26, 93, 186, 1209, 2418]>,
    <6, 455, [1, 6, 7, 10, 13, 15, 42, 70, 78, 91, 105, 130, 195, 546, 910, 1365]>,
    <210, 13, [1, 2, 21, 42, 65, 130, 1365, 2730]>,
    <3390, 1, [1, 5, 6, 30, 113, 565, 678, 3390]>,
    <714, 5, [1, 3, 10, 14, 17, 30, 35, 42, 51, 105, 170, 238, 510, 595, 714, 1785]>,
    <6, 649, [1, 2, 3, 6, 11, 22, 33, 59, 66, 118, 177, 354, 649, 1298, 1947, 3894]>,
    <4130, 1, [1, 2, 5, 7, 10, 14, 35, 59, 70, 118, 295, 413, 590, 826, 2065, 4130]>,
    <690, 7, [1, 6, 10, 14, 15, 21, 23, 35, 138, 210, 230, 322, 345, 483, 805, 4830]>
];

// Five further quotients that the filters leave undecided although an argument outside the code
// settles each of them as not hyperelliptic.  Each stays here until its argument is a stage.
//   * X_0^6(35)/W, two groups: the canonical quadrics of the stored models (0 relations) show that
//     the models are not hyperelliptic, but the models rest on the m = 0 term of Schofer's formula
//     at a composite level, so they are not on the list of proved models.
//   * X_0^14(51), X_0^38(15), X_0^158(3) with the groups below: at p = 3 the special fibre is the
//     genus-3 curve Z = X_0^D(N/3)/W' (W contains an involution swapping the two copies), Z is not
//     hyperelliptic by the dual graph at another prime, but Z mod 3 is not shown non-hyperelliptic
//     (the Weil polynomials of Z and its twists at 3 are all those of hyperelliptic curves).
open_gaps := [
    <6, 35, [1, 3, 5, 15], "model rests on the composite-level m = 0 term">,
    <6, 35, [1, 5, 7, 35], "model rests on the composite-level m = 0 term">,
    <14, 51, [1, 6, 14, 21, 34, 51, 119, 714], "special fibre at 3: Z = X_0^14(17)/{1,14,34,119} mod 3">,
    <38, 15, [1, 3, 10, 30, 38, 95, 114, 285], "special fibre at 3: Z = X_0^38(5)/{1,10,38,95} mod 3">,
    <158, 3, [1, 2, 237, 474], "special fibre at 3: Z = X_0^158(1)/{1,2} mod 3">
];
