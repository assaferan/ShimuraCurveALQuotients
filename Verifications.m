
intrinsic VerifyHHTable1(curves::SeqEnum[ShimuraQuot])
{Verify that Table 1 in [HH] (squarefree star curves) is reproduced when we reduce to a finite problem.
There are some discrepancies, see below.}
    Table1 := AssociativeArray([3..19]);
    Table1[3] := {97, 109, 113, 127, 128, 136, 139, 144, 149, 151,
          152, 162, 164, 169, 171, 175, 178, 179, 183, 185,
          187, 189, 194, 196, 203, 207, 217, 234, 236, 240,
          245, 246, 248, 249, 252, 258, 270, 282, 290, 294,
          295, 303, 310, 312, 315, 318, 329, 348, 420, 429,
          430, 455, 462, 476, 510};
    Table1[4] := {137, 148, 160, 172, 173, 176, 199, 200, 201, 202,
          214, 219, 224, 225, 228, 242, 247, 254, 259, 260,
          261, 262, 264, 267, 273, 275, 280, 300, 305, 306,
          308, 319, 321, 322, 335, 341, 342, 345, 350, 354,
          355, 366, 370, 374, 385, 399, 426, 434, 483, 546,
          570};
    Table1[5] := {157, 181, 192, 208, 212, 216, 218, 226, 235, 237,
          250, 253, 278, 279, 302, 323, 364, 371, 377, 378,
          391, 396, 402, 406, 410, 414, 418, 435, 438, 440,
          442, 444, 465, 494, 495, 595, 630, 714, 770, 798};
    Table1[6] := {163, 197, 211, 244, 265, 272, 274, 291, 297, 301,
          325, 336, 340, 470, 506, 561, 564, 690, 780, 858};
    Table1[7] := {193, 232, 268, 288, 296, 298, 309, 360, 372, 450,
          456, 460, 474, 492, 498, 504, 518, 558, 582, 660,
          870, 924};
    Table1[8] := {292, 408, 468, 480, 534, 540, 552, 606, 930, 966,
          990, 1020};
    Table1[9] := {516, 522, 528, 1110, 1140};
    Table1[10] := {600, 840, 1050, 1230, 1290};
    Table1[11] := {};
    Table1[12] := {2310};
    Table1[13] := {1260};
    Table1[14] := {2730};
    Table1[15] := {1470};
    Table1[16] := {};
    Table1[17] := {};
    Table1[18] := {};
    Table1[19] := {1680};
    // We add here a list of discrepancies
    // This is because the formula in [FH99] is slightly better than
    // the one in [HH96] accounting for all the cusps in the quotient
    Table1[3] diff:= {128};
    Table1[5] diff:= {495};
    Table1[6] diff:= {272, 297};
    Table1[7] diff:= {288, 296};
    // !! This is weird - N = 1170 should be appearing in Table 1,
    // but for some reason it does not.
    Table1[12] join:= {1170};
    assert Maximum([c`g : c in curves | c`D eq 1]) eq 19;
    for g in [3..19] do
    genus_g := {c`N : c in curves | (c`D eq 1) and (c`g eq g)};
    assert Table1[g] eq genus_g;
    end for;
end intrinsic;

intrinsic GetHHTable2() -> SeqEnum[RngIntElt]
{Returns Table 2 from [HH], with the exception of the genus 19 curve, which apparently can be ruled out by point counts!?}
    Table2 := [[] : g in [1..19]];
    Table2[3] := [ 127, 136, 144, 152, 162, 164, 171, 175, 183, 185,
           194, 196, 207, 217, 234, 240, 246, 252, 258, 270,
           282, 290, 294, 310, 312, 315, 318, 348, 420, 462,
           476, 510 ];
    Table2[4] := [ 160, 176, 264, 280, 300, 306, 322, 342, 345, 370,
           546, 570 ];
    Table2[5] := [ 216, 279, 396, 630, 714 ];
    Table2[6] := [ 336, 690 ];
    Table2[7] := [ 360, 450 ];
    Table2[10] := [840];
    Table2[12] := [2310];
    // This seems to be removed by considering p = 11, v = 2. Check!
    // Table2[19] := [1680];
    return Table2;
end intrinsic;

function GetModularByGenus(curves)
    g_max := Maximum([X`g : X in curves]);
    by_genus := [[X`N : X in curves | X`D eq 1 and
                      not assigned X`IsSubhyp and
                      X`g eq g] : g in [1..g_max]];
    return by_genus;
end function;

intrinsic VerifyHHTable2(starcurves::SeqEnum[ShimuraQuot])
{Verify that Table 2 in [HH] (squarefree star curves) is reproduced when we count points using trace formula on star curves.}
    Table2 := GetHHTable2();
    by_genus := GetModularByGenus(starcurves);
    assert Table2 eq by_genus;
end intrinsic;

intrinsic VerifyHHProposition1(starcurves::SeqEnum[ShimuraQuot])
{Verify that [HH, Prpoposition 1] rules out N = 194, 546 from the modular curve list.}
    Table2 := GetHHTable2();
    by_genus := GetModularByGenus(starcurves);
    Table2[3] := [N : N in Table2[3] | N ne 194];
    Table2[4] := [N : N in Table2[4] | N ne 546];
    assert Table2 eq by_genus;
end intrinsic;

intrinsic VerifyHasegawaTable3(starcurves::SeqEnum[ShimuraQuot])
{Verify that modular non-AL involution works for the curves in [Hasegawa, Table 3].
 !!! At the moment runs the check instead of testing the results of the filter. 
 !!! Change this when filter runs (should check if we lose something from FpAuts !?)}
    // Load the rows of the table, collected by value of N
    Table3 := [<136, ["V2"], 3, 0>, <144, ["V2", "V3", "V2 V3"], 3, 1>];
    Table3 cat:= [<152, ["V2"], 3, 1>, <171, ["V3"], 3, 0>, <207, ["V3"], 3, 2>, <234, ["V3"], 3, 1>];
    Table3 cat:= [<240, ["V2"], 3, 1>, <252, ["V3"], 3, 0>, <312, ["V2"], 3, 1>, <315, ["V3"], 3, 2>];
    Table3 cat:= [<160, ["V2"], 4, 2>, <176, ["V2"], 4, 0>, <264, ["V2"], 4, 1>, <280, ["V2"], 4, 1>];
    Table3 cat:= [<306, ["V3"], 4, 1>, <342, ["V3"], 4, 1>, <216, ["V2"], 5, 2>, <279, ["V3"], 5, 0>];
    Table3 cat:= [<396, ["V3"], 5, 3>, <630, ["V3"], 5, 2>, <336, ["V2"], 6, 2>, <360, ["V2"], 7, 2>, <450, ["V3"], 7, 3>];
    
    for row in Table3 do
        N, invs, g, gV := Explode(row);
        assert exists(X){X : X in starcurves | X`D eq 1 and X`N eq N};
        is_hyp, inv, fix_v := CheckModularNonALInvolutionModSym(X);
        assert X`g eq g;
        if is_hyp ne -1 then
            invs2 := [v cat " W1" : v in invs];
            assert assigned inv;
            if (inv notin (invs cat invs2)) then
                vprintf ShimuraQuotients, 3: "row = ", row, "inv = ", inv, ", invs = ", invs;
                break;
            end if;
            assert (inv in invs) or (inv in invs2);
        end if;
        if is_hyp eq 0 then
            assert assigned fix_v;
            assert fix_v eq 2*g - 4*gV + 2;
        elif is_hyp eq 1 then
            assert gV eq 0;
        end if;
    end for;

    return;
end intrinsic;

intrinsic VerifyFHTheorem3(curves::SeqEnum)
    {Verify Theorem 3 of [FH99]: there are exactly 32 pairs (N, W') for which
    X_0(N)/W' is a hyperelliptic curve of genus g >= 3 whose hyperelliptic
    involution v is of Atkin-Lehner type. Run after UpdateByGenus on the full
    quotient list.}
    // Entries are <N, generators of W', g, v>, transcribed from the table
    // on p. 110 of [FH99].
    Table3 := [
        <46, {2}, 3, 23>,
        <51, {3}, 3, 17>,
        <55, {5}, 3, 11>,
        <56, {8}, 3, 7>,
        <60, {4}, 3, 15>,
        <60, {12}, 4, 20>,
        <60, {60}, 3, 4>,
        <62, {2}, 4, 31>,
        <66, {6}, 4, 11>,
        <66, {66}, 3, 6>,
        <69, {3}, 4, 23>,
        <70, {10}, 4, 14>,
        <70, {14}, 3, 10>,
        <78, {6}, 6, 26>,
        <78, {26}, 3, 6>,
        <78, {2, 3}, 3, 13>,
        <78, {13, 6}, 3, 3>,
        <87, {3}, 5, 29>,
        <92, {4}, 5, 23>,
        <92, {92}, 4, 4>,
        <94, {2}, 6, 47>,
        <94, {94}, 4, 2>,
        <95, {5}, 5, 19>,
        <95, {19}, 3, 5>,
        <105, {3, 5}, 3, 7>,
        <105, {3, 7}, 3, 5>,
        <105, {7, 15}, 3, 3>,
        <110, {2, 5}, 4, 11>,
        <110, {2, 11}, 3, 5>,
        <110, {5, 22}, 3, 2>,
        <119, {7}, 6, 17>,
        <119, {17}, 4, 7>
    ];
    expected := {<t[1], AllALsFromGens(t[2], t[1]), t[3],
                  AllALsFromGens(t[2] join {t[4]}, t[1])> : t in Table3};
    assert #expected eq 32;
    // X_0(N)/W' has an AL involution as hyperelliptic involution iff
    // some quotient by an index-2 overgroup <W', v> has genus 0
    actual := {};
    for c in curves do
        if (c`D ne 1) or (c`g lt 3) then continue; end if;
        for j in c`Covers do
            if curves[j]`g eq 0 then
                // UpdateByGenus must have marked it hyperelliptic
                assert c`IsHyp;
                Include(~actual, <c`N, c`W, c`g, curves[j]`W>);
            end if;
        end for;
    end for;
    // Theorem 3 only covers nontrivial W'; the curves with W' = {1} are
    // the hyperelliptic X_0(N) themselves ([FH99] Theorem 1, due to Ogg).
    // Of those, exactly the following have g >= 3 with hyperelliptic
    // involution of AL type (in particular X_0(37), X_0(40), X_0(48),
    // whose hyperelliptic involutions are not AL, do not show up here).
    expected1 := {<30, 3, {1, 15}>, <33, 3, {1, 11}>, <35, 3, {1, 35}>,
                  <39, 3, {1, 39}>, <41, 3, {1, 41}>, <46, 5, {1, 23}>,
                  <47, 4, {1, 47}>, <59, 5, {1, 59}>, <71, 6, {1, 71}>};
    assert expected1 eq {<t[1], t[3], t[4]> : t in actual | t[2] eq {1}};
    assert expected eq {t : t in actual | t[2] ne {1}};
end intrinsic;

intrinsic VerifyFHProposition4(curves::SeqEnum)
    {Verify that FindIsomorphicCurveProp4 finds the isomorphisms listed
    after Proposition 4 of [FH99] (p. 114). The list there is not
    comprehensive, so we only check that each listed isomorphism is found.}
    // Entries are <M, generators of W on X_0(M), N = 2M, generators of W'
    // on X_0(N)> with X_0(M)/W isomorphic to X_0(N)/W'. The curves listed
    // as X_0^*(N) in the paper are spelled out by generators of W(N).
    isos := [
        <34, {Integers()|}, 68, {4}>,
        <82, {41}, 164, {4, 41}>,
        <42, {7}, 84, {4, 7}>,
        <98, {49}, 196, {4, 49}>,
        <106, {53}, 212, {4, 53}>,
        <118, {59}, 236, {4, 59}>,
        <154, {7, 11}, 308, {4, 7, 11}>,
        <174, {3, 29}, 348, {4, 3, 29}>,
        <90, {9}, 180, {4, 9}>,
        <90, {5}, 180, {4, 5}>,
        <90, {45}, 180, {4, 45}>,
        <198, {9, 11}, 396, {4, 9, 11}>,
        <102, {51}, 204, {4, 51}>,
        <210, {3, 5, 7}, 420, {4, 3, 5, 7}>,
        <238, {7, 17}, 476, {4, 7, 17}>,
        <138, {23}, 276, {4, 23}>
    ];
    lut := AssociativeArray();
    for i in [1..#curves] do
        lut[<curves[i]`D, curves[i]`N, curves[i]`W>] := i;
    end for;
    for iso in isos do
        M, gensM, N, gensN := Explode(iso);
        X := curves[lut[<1, M, AllALsFromGens(gensM, M)>]];
        Y := FindIsomorphicCurveProp4(X, curves, lut);
        assert Type(Y) eq ShimuraQuot;
        assert Y`N eq N;
        assert Y`W eq AllALsFromGens(gensN, N);
        // isomorphic curves have the same genus
        assert X`g eq Y`g;
    end for;
end intrinsic;

intrinsic VerifyFHProposition5(curves::SeqEnum)
    {Verify that FindIsomorphicCurveProp5 finds the isomorphisms listed
    after Proposition 5 of [FH99] (p. 114). The list there is not
    comprehensive, so we only check that each listed isomorphism is found.}
    // Entries are <N, generators of W', generators of W''> with
    // X_0(N)/W' isomorphic to X_0(N)/W''.
    isos := [
        <90, {18, 10}, {2, 5}>,
        <99, {99}, {11}>,
        <180, {36, 20}, {5, 36}>,
        <198, {11, 18}, {2, 99}>
    ];
    lut := AssociativeArray();
    for i in [1..#curves] do
        lut[<curves[i]`D, curves[i]`N, curves[i]`W>] := i;
    end for;
    for iso in isos do
        N, gensW, gensW2 := Explode(iso);
        X := curves[lut[<1, N, AllALsFromGens(gensW, N)>]];
        Y := FindIsomorphicCurveProp5(X, curves, lut);
        assert Type(Y) eq ShimuraQuot;
        assert Y`N eq N;
        assert Y`W eq AllALsFromGens(gensW2, N);
        // isomorphic curves have the same genus
        assert X`g eq Y`g;
    end for;
end intrinsic;

intrinsic VerifyFHTable3(curves::SeqEnum)
    {Verify that FilterByComplicatedALFixedPointsOnQuotient (implementing
    [FH99] Prop. 6) proved non-hyperellipticity for all 34 pairs (N, W')
    of [FH99] Table 3 (p. 116), by checking the actual curve data.}
    // Entries are <N, generators of W'>
    Table3 := [
        <58, {2}>, <58, {58}>,
        <76, {19}>, <76, {76}>,
        <86, {86}>,
        <102, {2, 17}>, <102, {17, 6}>,
        <106, {106}>,
        <114, {2, 57}>,
        <122, {122}>,
        <124, {31}>,
        <130, {2, 65}>,
        <132, {3, 44}>, <132, {11, 12}>,
        <134, {134}>,
        <140, {7, 20}>, <140, {20, 28}>,
        <150, {2, 75}>, <150, {3, 50}>,
        <170, {2, 17}>, <170, {10, 34}>,
        <174, {2, 87}>,
        <182, {14, 26}>,
        <186, {6, 62}>,
        <190, {2, 95}>, <190, {10, 38}>,
        <198, {2, 99}>,
        <204, {3, 68}>, <204, {12, 68}>,
        <210, {5, 6, 14}>,
        <222, {6, 74}>,
        <230, {5, 46}>,
        <330, {5, 11, 6}>,
        <390, {6, 10, 26}>
    ];
    assert #Table3 eq 34;
    lut := AssociativeArray();
    for c in curves do
        lut[<c`D, c`N, c`W>] := c;
    end for;
    for t in Table3 do
        N := t[1]; W := AllALsFromGens(t[2], N);
        if not IsDefined(lut, <1, N, W>) then
            error Sprintf("Table 3 curve D=1, N=%o, W=%o not found in curves", N, W);
        end if;
        c := lut[<1, N, W>];
        if not (assigned c`IsHyp and not c`IsHyp) then
            error Sprintf("Table 3 curve D=1, N=%o, W=%o not proved non-hyperelliptic", N, W);
        end if;
    end for;
end intrinsic;