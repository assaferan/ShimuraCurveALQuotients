
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
                print "row = ", row, "v = ", v, ", invs = ", invs;
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