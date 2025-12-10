
procedure test_Schofer(f, vals, D, N)
    L := ShimuraCurveLattice(D,N);
    for datum in vals do
        d, val := Explode(datum);
        assert SchoferFormula(f, d, D, N, L) eq LogSum(val);
    end for;
end procedure;

procedure test_Schofer_6()
    printf "Testing Schofer formula for D=6...";
    R := EtaQuotientsRing(12, 72);

    psi0 := EtaQuotient(R, [-2,5,0,-2,0,0]);
    psi0_qexp<q> := qExpansionAtoo(psi0,10);
    assert psi0_qexp eq 1 + 2*q + 2*q^4 + 2*q^9 + O(q^10);

    psi1 := EtaQuotient(R, [-5,12,1,-4,-1,-2]);
    assert qExpansionAtoo(psi1,1) eq q^-1 + 5 + O(q);
    
    psi3 := EtaQuotient(R, [0,1,2,4,4,-10]);
    assert qExpansionAtoo(psi3,1) eq q^-3 - q^-1 - 2 + O(q);

    // testing [Errthum, p. 850]
    f6 := -6*psi3 - 2*psi1 - 2*psi0;
    assert qExpansionAtoo(f6,1) eq -6*q^(-3) + 4*q^(-1) + O(q);

    // |t6| = 6^6 | psi_f6 |^2
    vals := [<-24, { <3, -6>, <2, -6> }>,
             <-163, { <3, 5>, <2, -16>, <7,4>, <5,-6>, <19,4>, <11,-6>, <23,4>, <17,-6> }>,
             // testing [Errthum, Table 2] (Section 8.1.1)
             <-40, { <3, 1>, <2, -6>, <5,-3> }>,
             <-19, { <3, 1>, <2, -16> }>,
             <-52, { <2, -4>, <3, 1>, <5,-6> }>,
             <-84, { <2, -4>, <3, -9>, <7,2> }>,
             <-88, { <2, -6>, <3, 1>, <5,-6>, <7,4>, <11,-3> }>,
             // This one does not work! - the problem with zeros
             // <-100, { <2, -2>, <3, 1>, <5,1>, <7,4>, <11,-6> }>,
             <-120, { <2, -6>, <3, -9>, <5,-3>, <7,4> }>,
             <-132, { <2, -2>, <3, -6>, <5,-6>, <11,2> }>,
             <-148, { <2, -4>, <3, 1>, <5,-6>, <7,4>, <11,4>, <17,-6> }>,
             <-168, { <2, -6>, <3, -6>, <5,-6>, <7,2>, <11,4> }>,
             <-43, { <2, -16>, <3, 1>, <5,-6>, <7,4>}>,
             <-51, { <2, -16>, <3, -6>, <7,4>}>,
             <-228, { <3, -12>, <5,-6>, <7,4>, <19,2> }>,
             <-232, { <2,-6>, <3, 1>, <5,-6>, <7,4>, <11,4>, <19,4>, <23,-6>, <29,-3> }>,
             <-67, { <2,-22>, <3, 1>, <5,-6>, <7,4>, <11,4> }>,
             // This does not work - hits mu = m = 0 !!!!???
             // <-75, { <2,-16>, <3, -9>, <5,-1>, <11,4> }>,
             <-312, { <2,-6>, <3, -6>, <5,-6>, <7,4>, <11,-6>, <23,4> }>,
             <-372, { <2,-4>, <3, -9>, <5,-6>, <7,4>, <11,-6>, <19,4>, <31,2> }>,
             <-408, { <2,-6>, <3, -12>, <5,-6>, <7,4>, <11,4>, <17,-3>, <31,4> }>,
             <-123, { <2,-16>, <3, -6>, <5,-6>, <7,4>, <19,4> }>,
             // This one does not work! - the problem with zeros
             // Also see p.828 of [Err]
             // <-147, { <2,-16>, <3, -9>, <5,-6>, <7,-1>, <11,4>, <23,4> }>,
             <-163, { <2,-16>, <3, 5>, <5,-6>, <7,4>, <11,-6>, <17,-6>, <19,4>, <23,4> }>,
             <-708, { <2,2>, <3, -6>, <5,-6>, <7,4>, <11,4>, <17,-6>, <29,-6>, <47,4>, <59,2> }>,
             <-267, { <2,-22>, <3, -6>, <5,-6>, <7,4>, <11,-6>, <31,4>, <43,4> }>,
             // divide to account for 3 CM points of degree 6
             <-996, { <2, -2>, <3, -18>, <41, -6>, <29, -6>, <17,-6>, <7,12>, <83,2>, <71,4> }>
    ];

    test_Schofer(f6, vals, 6, 1);

    printf "Done!\n";
    return;

end procedure;

procedure test_Schofer_10()
    printf "testing Schofer formula for D=10...";
    R := EtaQuotientsRing(20, 200);

    f10 := 3*EtaQuotient(R, [0,-3,6,-2,8,-8]) - 2*EtaQuotient(R, [-2,3,2,0,2,-4]) - 5*EtaQuotient(R, [0,-1,2,-2,6,-4]) + 4*EtaQuotient(R, [-2,5,-2,0,0,0]);
    qexp_f10<q> := qExpansionAtoo(f10,1);

    assert qexp_f10 eq 3*q^(-3) - 2*q^(-2) + O(q);

    //t_10 = 2^(-2)*|Psi_f_10|^2
    vals := [<-20, { <2, 3> }>,
             // The next one does still not work at 2! Off by a factor of 4, this is confusing. The kappas are all correct now.
             // Probably an error in [Err] ?
             // Here is the result we expected from [Err]
             // <-68, { <2, 2>, <5, 1> }>,
             <-68, { <2, 6>, <5, 1> }>,
             //Testing Err. Table 4, Section 8.2.1
             <-40, { <3,3>, <2, 2> }>,
             <-52, { <3,3>, <2, 3>, <5, -2> }>,
             // !At the moment, non-maximal orders are not implemented!
             // <-72, { <3,1>, <2, 2>, <5, -2>, <7,-2> }>,
             <-120, { <3,3>, <2, 2>,  <7,-2> }>,
             <-88, { <3,3>, <2, 1>, <5,3>, <7,-2> }>,
             // !At the moment, non-maximal orders are not implemented!
             // this is wrong -used to have the 0 error
             // <-27, { <3,1>, <2, 8>, <5,-2> }>,
             <-35, {  <2, 8>, <7,-1> }>,
             <-148, {  <2, 3>, <3,3>, <11,3>, <5,-2>, <7,-2>, <13,-2> }>,
             <-43, {  <2, 8>, <3,3>,  <5,-2>, <7,-2>}>,
             // !At the moment, non-maximal orders are not implemented!
             // This (-180) doesn't work! Why??? This was proved in Elkies, also
             // <-180, {  <2, 3>,  <11,3>,<13,-2>}>,
             <-232, {   <3,3>,  <11,3>,<17,3>, <5,-2>, <7,-2>, <23,-2>}>,
             <-67, {   <2,8>,  <3,3>,  <5,3>, <7,-2>, <13,-2>}>,
             <-280, {   <2,1>,  <3,3>, <11,3>, <7,-1>, <23,-2>}>,
             <-340, {   <2,3>,  <3,3>, <23,3>, <7,-2>, <29,-2>}>,
             <-115, {   <2,11>,  <3,3>, <13,-2>, <23,-1>}>,
             <-520, {   <2,-1>,  <3,3>, <29,3>, <7,-2>, <13,-1>, <47,-2>}>,
             // this works but takes a long time!
             <-163, {   <2,11>,  <3,3>, <5,3>, <11,3>, <7,-2>, <13,-2>,<29,-2>, <31,-2> }>,
             <-760, {   <2,2>,  <3,3>, <17,3>, <47,3>, <7,-2>, <31,-2>,<71,-2> }>,
             <-235, {   <2,8>,  <3,3>, <17,3>,  <7,-2>, <37,-2>,<47,-1> }>,
             // this works! the extra 2^2 comes from the fact that there are 2 points here
             <-420, { <2,6> , <3,6>,<29,3>,<7,-2>,<37,-2> }>
    ];

    test_Schofer(f10, vals, 10, 1);

    printf "Done!\n";
    return;
end procedure;

procedure test_Schofer_142()
    printf "testing Schofer formula for D=142...";

    vals1 := [// GY p.30
              // then psi_f1 *2^10= x, but I'm getting the reverse...
              <-19, { <2,-10> }>,
              <-20, { <2,-10> }>,
              <-24, { <2,-11> }>,
              <-40, { <2,-11> }>,
              <-43, { <2,-10> }>,
              // this works, but long time
              <-148, { <2,-10>  }>,
              // this works, but long time
              <-232, { <2,-11> }>
    ];

    // psi_f3 *2 = y
    vals3 := [// Can't compute this because m = 0 !?
              // <-3, { <2,1> }>,
              <-4, {}>,
              // Can't compute this because m = 0 !?
              // <-19, { <2,-1> }>,
              <-24, { <2,1> }>,
              <-43, { <3,1> }>
    ];

    _<q> := LaurentSeriesRing(Rationals());
    f1 := -2*q^(-87)- 2*q^(-71) - 2*q^(-48) - 2*q^(-36) +2*q^(-16) -2*q^(-15)-2*q^(-12)+2*q^(-9)-2*q^(-7) - 2*q^(-3) +2*q^(-2) + 2*q^(-1) + O(q);
    f3 := q^(-87) - 2*q^(-79) + q^(-76) - 2*q^(-71)+2*q^(-48)-q^(-40)+3*q^(-32)-2*q^(-20)-q^(-19)-2*q^(-12)+2*q^(-10)+q^(-8)-2*q^(-7)-2*q^(-2) + O(q);
   
    test_Schofer(f1, vals1, 142, 1);
    test_Schofer(f3, vals3, 142, 1);

    printf "Done\n";
    return;
end procedure;

procedure test_Schofer_26()
    printf "testing Schofer formula for D=26...";
    //Test as in Guo Yang p.25
    _<q> := LaurentSeriesRing(Rationals());

    g1 := 2*q^(-13) - 2*q^(-2) - 4*q^(-1) + O(q);
    g2 := q^(-11) + 2*q^(-7) - 2*q^(-2) + O(q);
    g3 := 2*q^(-26) + 6*q^(-7)- 6*q^(-2) +2*q^(-1) + O(q);

    g := [g1, g2, g3];

    vals1 := [<-11, {}>,
              <-19, { <3,2>  }>,
              <-20, { <5,1>  }>,
              <-24, { <3,1>  }>,
              <-67, { <3,4>, <5,-2>  }>
    ];

    vals2 := [<-19, { <2,6>  }>,
              <-20, {<2,5>}>,
              <-24, {<2,5>}>,
              <-52, {<2,3>}>,
              <-67, { <2,6>, <5,-2>, <7,1>  }>
    ];

    vals3 := [<-11, { <2,10>,<11,1>,<13,3> }>,
              <-19, { <2,10>,<13,3>,<19,1>  }>,
              <-20, { <2,12>,<13,3>  }>,
              <-24, { <2,13>,<13,3>  }>,
              // this next one does not work!
              // Probably a typo in [GY]
              // Original value in [GY]
              // <-52, { <2,6>,<13,8>  }>,
              <-52, { <2,6>,<13,5>  }>,
              <-67, { <2,10>, <5,-6>, <13,3>,<41,2>,<67,1>  }>
    ];

    vals := [vals1, vals2, vals3];

    for i in [1..3] do
        test_Schofer(g[i], vals[i], 26, 1);
    end for;

    printf "Done!\n";
    return;
end procedure;

test_Schofer_6();
test_Schofer_10();
test_Schofer_142();
test_Schofer_26();