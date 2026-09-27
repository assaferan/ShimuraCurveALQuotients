// tests/CMFieldsOfDefinition.m
//
// CI checks for CMFieldsOfDefinition.m and the scan RationalandQuadraticCMPoints that uses it.
//
//   [1] CMClassLists.  SOURCE: [Kl] J. Klaise, "Orders in quadratic imaginary fields of small
//       class number", Warwick report (2013), table p. 19, and [Wat] M. Watkins, "Class numbers
//       of imaginary quadratic fields", Math. Comp. 73 (2004), table p. 936.  A brute force over
//       Magma's ClassNumber must reproduce both; CMClassLists must then equal it.
//   [2] DegreeOfFieldOfDefinitionOfCMPoint vs the degrees of FieldsOfDefinitionOfCMPointFast on
//       every W of three bases.  NO EXTERNAL SOURCE.  The two share the existence test, D_R, N_R,
//       N*_R, W_gal and cc_active, so this checks only the independent parts: #W_gal / 2^eps
//       against the Artin-class subgroup #alSub, h(R) against the norm group, and the cc halving
//       against the constructed fixed field.  The degree is checked against [GY] in tests/CMPoints.m.
//   [3] RationalCMDiscs and RationalandQuadraticCMPoints on X_0^6(1)* and X_0^10(1)*.
//       SOURCE: [Err] E. Errthum, "Singular moduli of Shimura curves", Canad. J. Math. 63 (2011),
//       arXiv:0711.4316v2, Table 2 (D = 6, 27 rows), Table 4 (D = 10, 22 rows), Sec. 7.4 (the
//       irrational point -68 on X_10^*).  Asserted: (i) every table row is rational; (ii) the
//       fundamental rational discriminants are exactly the table's (Prop. 2.4 covers maximal
//       orders only); (iii) for D = 6 only, full completeness with one point each (Sec. 2.4: "17
//       of the 27 rational CM points").  [Err] gives no total for D = 10 and no degrees.

// [1]
procedure test_ClassNumberTables()
    printf "Testing CMClassLists against [Kl] and brute force...";
    // [Kl, p. 19]: <h, number of orders, largest |disc|>
    Kl := [<1, 13, 163>, <2, 29, 427>, <3, 25, 907>, <4, 84, 1555>, <5, 29, 2683>,
           <6, 101, 4075>, <7, 38, 5923>, <8, 208, 7987>, <16, 531, 35275>];
    // [Wat, p. 936]: <h, number of fundamental discriminants, largest |disc|>
    Wat := [<1, 9, 163>, <2, 18, 427>, <4, 54, 1555>, <8, 131, 6307>, <16, 322, 31243>];
    bound := 40000;
    assert bound gt Max([t[3] : t in Kl]);
    all_orders := AssociativeArray();
    fundamental := AssociativeArray();
    for h in [1..16] do
        all_orders[h] := {Integers() | };
        fundamental[h] := {Integers() | };
    end for;
    for n in [3..bound] do
        d := -n;
        if d mod 4 notin {0, 1} then continue; end if;
        h := ClassNumber(d);
        if h gt 16 then continue; end if;
        Include(~all_orders[h], d);
        if IsFundamentalDiscriminant(d) then Include(~fundamental[h], d); end if;
    end for;
    for t in Kl do
        h, count, largest := Explode(t);
        assert #all_orders[h] eq count and Max([-d : d in all_orders[h]]) eq largest;
    end for;
    for t in Wat do
        h, count, largest := Explode(t);
        assert #fundamental[h] eq count and Max([-d : d in fundamental[h]]) eq largest;
    end for;

    CNs := CMClassLists();
    assert Keys(CNs) eq {1..8};
    for h in [1..8] do
        assert CNs[h] eq all_orders[h];
    end for;
    printf "Done!\n";
end procedure;

// [2]
procedure test_DegreeMatchesFields()
    CNs := CMClassLists();
    discs := Sort(Setseq(&join[CNs[h] : h in [1..4]]));
    for DN in [<6, 1>, <1, 30>, <6, 5>] do
        D, N := Explode(DN);
        printf "Cross-checking the CM degree against the fields on every W for (D,N) = (%o,%o)...", D, N;
        nchecks := 0;
        for Wdata in ALSubgroups(D*N) do
            X := CreateShimuraQuot(D, N, Wdata[1]);
            for d in discs do
                deg := DegreeOfFieldOfDefinitionOfCMPoint(X, d);
                Fs := FieldsOfDefinitionOfCMPointFast(X, d);
                if deg eq 0 then
                    assert #Fs eq 0;
                else
                    assert #Fs gt 0;
                    assert &and[AbsoluteDegree(F) eq deg : F in Fs];
                end if;
                nchecks +:= 1;
            end for;
        end for;
        // guard against the loop silently checking nothing
        assert nchecks eq #ALSubgroups(D*N) * #discs;
        printf "Done! (%o checks)\n", nchecks;
    end for;
end procedure;

// [3]
procedure test_ErrthumTables()
    // [Err, arXiv:0711.4316v2, Table 2, Sec. 8.2, "Coordinates of Rational CM Points on X_6^*"],
    // 27 rows, in the table's row order.
    rat6 := [-3, -4, -24, -40, -52, -19, -84, -88, -100, -120, -132, -148, -168, -43, -51, -228,
             -232, -67, -75, -312, -372, -408, -123, -147, -163, -708, -267];
    // [Err, arXiv:0711.4316v2, Table 4, Sec. 8.5, "Coordinates of Rational CM Points on
    // X_10^*"], 22 rows, in the table's row order.
    rat10 := [-3, -8, -20, -40, -52, -72, -120, -88, -27, -35, -148, -43, -180, -232, -67, -280,
              -340, -115, -520, -163, -760, -235];
    // guard the transcription: row counts as printed, no repeated rows
    assert #rat6 eq 27 and #Set(rat6) eq 27;
    assert #rat10 eq 22 and #Set(rat10) eq 22;
    // [Err, Sec. 7.4]: "the irrational CM point with discriminant -68" on X_10^*
    irr := AssociativeArray();
    irr[6] := {Integers() | };
    irr[10] := {-68};
    // [Err, Sec. 2.4]: "Elkies was able to compute the coordinates of 17 of the 27 rational CM
    // points (see Table 2)."  Only D = 6 has a stated total.
    total_known := AssociativeArray();
    total_known[6] := true;
    total_known[10] := false;
    fund := func<S | {d : d in S | IsFundamentalDiscriminant(d)}>;
    for datum in [<6, Set(rat6)>, <10, Set(rat10)>] do
        D, rat := Explode(datum);
        printf "Testing the rational CM points on X_0^%o(1)* against [Err]...", D;
        X := CreateShimuraQuot(D, 1, Set(Divisors(D)));
        rat_tab := RationalCMDiscs(X);
        found := Keys(rat_tab);
        // (i) membership: every table row is a rational CM point
        assert rat subset found;
        // (ii) completeness at maximal orders.  [Err, Prop. 2.4]: "P_Delta is a rational point on
        // X_D^* if and only if the class group of k is generated by ideals I subset R such that
        // I^2 = (p) for some p | D", with R "the maximal order in the quadratic imaginary field of
        // discriminant Delta"; then "In the case of d(B) = 2, all such fields are known, and thus
        // the rational CM points can be identified. (See Table 2 for D = 6 and Table 4 for
        // D = 10.)"
        assert fund(found) eq fund(rat);
        // (iii) full completeness, one point per discriminant: only where [Err] states the total
        // (27 points, 27 distinct Delta in Table 2).
        if total_known[D] then
            assert found eq rat;
            assert &and[rat_tab[d] eq 1 : d in rat];
        end if;
        for d in irr[D] do
            // a CM point exists (degree > 0) and it is not rational
            assert DegreeOfFieldOfDefinitionOfCMPoint(X, d) gt 1;
            assert d notin found;
        end for;
        // The scan's table holds only fundamental discriminants, so it must find every fundamental
        // one of [Err], and none of them among the degree-2 points.
        rat_scan, quad_scan := RationalandQuadraticCMPoints(X : bd := 8, coprime_to_level := false);
        rat_scan := {p[1] : p in rat_scan};
        quad_scan := {p[1] : p in quad_scan};
        assert fund(rat) subset rat_scan and #(fund(rat) meet quad_scan) eq 0;
        assert #(irr[D] meet rat_scan) eq 0;
        printf "Done! (%o rational, %o fundamental)\n", #rat, #fund(rat);
    end for;
end procedure;

test_ClassNumberTables();
test_DegreeMatchesFields();
test_ErrthumTables();
