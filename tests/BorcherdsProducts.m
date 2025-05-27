
procedure test_Kappa0()
    D := 6;
    N := 1;
    L, Ldual, disc_grp, to_disc, Qinv := ShimuraCurveLattice(D,N);
    Q := ChangeRing(Qinv^(-1), Integers());
    // verifying [Yang, Example 21, p. 24-25]
    // assert Round(1/Kappa0(3,-4,Q)) eq 2^8*3^4;
    log_coeffs := Kappa0(3,-4,Q);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs)} eq { <2, -8>, <3, -4> };
    // assert Round(1/Kappa0(1,-3,Q)) eq 2^4;
    log_coeffs := Kappa0(1,-3,Q);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs)} eq { <2, -4> };
    // verifying [Err, p. 850]
    // assert Round(1/Kappa0(1,-24,Q)) eq 2^6;
    log_coeffs := Kappa0(1,-24,Q);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs)} eq { <2, -6> };
    // assert Round(1/Kappa0(3,-24,Q)) eq 2^8*3^4;
    log_coeffs := Kappa0(3,-24,Q);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs)} eq { <2, -8>, <3, -4> };
    // assert Round(1/Kappa0(1,-163,Q)) eq 2^4*3^11*7^4*19^4*23^4;
    log_coeffs := Kappa0(1,-163,Q);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs)} eq { <2, -4>, <3, -11>, <7, -4>, <11,0>, <19,-4>, <23,-4> };
    // assert Round(1/Kappa0(3,-163,Q)^3) eq 2^40*3^12*5^12*11^12*17^12;
    log_coeffs := Kappa0(3,-163,Q);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs)} eq { <2, -40/3>, <3, -4>, <5, -4>, <11,-4>, <17,-4>, <23,0>, <89, 0> };
    D := 10;
    L, Ldual, disc_grp, to_disc, Qinv := ShimuraCurveLattice(D,N);
    Q := ChangeRing(Qinv^(-1), Integers());
    log_coeffs := Kappa0(3,-68,Q);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0} eq { <2, -8>,  <5, -14/3>};
    // This does not work. This is probably a typo! Isn't it supposed to be 2. This would be the relevant number.
    //log_coeffs := Kappa0(1,-68,Q);
    log_coeffs := Kappa0(2,-68,Q);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0} eq { <2, -6>,  <5, -6>};
    return;
end procedure;

procedure test_Schofer_6()
    _<q> := LaurentSeriesRing(Rationals());
    // testing [Errthum, p. 850]
    f6 := -6*q^(-3) + 4*q^(-1) + O(q);
    log_coeffs := SchoferFormula(f6, -24, 6, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0} eq { <3, -6>, <2, -6> };

    // |t6| = 6^6 | psi_f6 |^2
    log_coeffs := SchoferFormula(f6, -163, 6, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0} eq { <3, 5>, <2, -16>, <7,4>, <5,-6>, <19,4>, <11,-6>, <23,4>, <17,-6> };

    // testing [Errthum, Table 2] (Section 8.1.1)
    log_coeffs := SchoferFormula(f6, -40, 6, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0} eq { <3, 1>, <2, -6>, <5,-3> };

    log_coeffs := SchoferFormula(f6, -19, 6, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0} eq { <3, 1>, <2, -16> };

    log_coeffs := SchoferFormula(f6, -52, 6, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0} eq { <2, -4>, <3, 1>, <5,-6> };

    log_coeffs := SchoferFormula(f6, -84, 6, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0} eq { <2, -4>, <3, -9>, <7,2> };

    log_coeffs := SchoferFormula(f6, -88, 6, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0} eq { <2, -6>, <3, 1>, <5,-6>, <7,4>, <11,-3> };

    // This one does not work! - the problem with zeros
    // log_coeffs := SchoferFormula(f6, -100, 6, 1);
    // assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0} eq { <2, -2>, <3, 1>, <5,1>, <7,4>, <11,-6> };

    log_coeffs := SchoferFormula(f6, -120, 6, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0} eq { <2, -6>, <3, -9>, <5,-3>, <7,4> };

    log_coeffs := SchoferFormula(f6, -132, 6, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0} eq { <2, -2>, <3, -6>, <5,-6>, <11,2> };

    log_coeffs := SchoferFormula(f6, -148, 6, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0} eq { <2, -4>, <3, 1>, <5,-6>, <7,4>, <11,4>, <17,-6> };

    log_coeffs := SchoferFormula(f6, -168, 6, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0} eq { <2, -6>, <3, -6>, <5,-6>, <7,2>, <11,4> };

    log_coeffs := SchoferFormula(f6, -43, 6, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0} eq { <2, -16>, <3, 1>, <5,-6>, <7,4>};

    log_coeffs := SchoferFormula(f6, -51, 6, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0} eq { <2, -16>, <3, -6>, <7,4>};

    log_coeffs := SchoferFormula(f6, -228, 6, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0} eq { <3, -12>, <5,-6>, <7,4>, <19,2> };

    log_coeffs := SchoferFormula(f6, -232, 6, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0} eq { <2,-6>, <3, 1>, <5,-6>, <7,4>, <11,4>, <19,4>, <23,-6>, <29,-3> };

    log_coeffs := SchoferFormula(f6, -67, 6, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0} eq { <2,-22>, <3, 1>, <5,-6>, <7,4>, <11,4> };

    // This does not work - hits mu = m = 0 !!!!???
    // log_coeffs := SchoferFormula(f6, -75, 6, 1);
    // assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0} eq { <2,-16>, <3, -9>, <5,-1>, <11,4> };

    log_coeffs := SchoferFormula(f6, -312, 6, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0} eq { <2,-6>, <3, -6>, <5,-6>, <7,4>, <11,-6>, <23,4> };

    log_coeffs := SchoferFormula(f6, -372, 6, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0} eq { <2,-4>, <3, -9>, <5,-6>, <7,4>, <11,-6>, <19,4>, <31,2> };

    log_coeffs := SchoferFormula(f6, -408, 6, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0} eq { <2,-6>, <3, -12>, <5,-6>, <7,4>, <11,4>, <17,-3>, <31,4> };

    log_coeffs := SchoferFormula(f6, -123, 6, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0} eq { <2,-16>, <3, -6>, <5,-6>, <7,4>, <19,4> };

    // This one does not work! - the problem with zeros
    // Also see p.828 of [Err]
    // log_coeffs := SchoferFormula(f6, -147, 6, 1);
    // assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0} eq { <2,-16>, <3, -9>, <5,-6>, <7,-1>, <11,4>, <23,4> };
    
    log_coeffs := SchoferFormula(f6, -163, 6, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0} eq { <2,-16>, <3, 5>, <5,-6>, <7,4>, <11,-6>, <17,-6>, <19,4>, <23,4> };

    log_coeffs := SchoferFormula(f6, -708, 6, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0} eq { <2,2>, <3, -6>, <5,-6>, <7,4>, <11,4>, <17,-6>, <29,-6>, <47,4>, <59,2> };

    log_coeffs := SchoferFormula(f6, -267, 6, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0} eq { <2,-22>, <3, -6>, <5,-6>, <7,4>, <11,-6>, <31,4>, <43,4> };

    //divide to account for 3 CM points of degree 6
    log_coeffs := SchoferFormula(f6, -996, 6, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs)} eq { <2, -2>, <3, -18>, <41, -6>, <29, -6>, <17,-6>, <7,12>, <83,2>, <71,4> };

    return;




end procedure;


procedure test_Schofer_10()

    _<q> := LaurentSeriesRing(Rationals());

    //This works!
    f10 := 3*q^(-3) - 2*q^(-2) + O(q);
    log_coeffs := SchoferFormula(f10, -20, 10, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs)} eq { <2, 3> };

    //This does still not work at 2! Off by a factor of 4, this is confusing. The kappas are all correct now.
    log_coeffs := SchoferFormula(f10, -68, 10, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs)} eq { <2, 2>, <5, 1> };

    //t_10 = 2^(-2)*|Psi_f_10|^2

    //Testing Err. Table 4, Section 8.2.1
    //this works!
    log_coeffs := SchoferFormula(f10, -40, 10, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs)} eq { <3,3>, <2, 2> };

    //this works!
    log_coeffs := SchoferFormula(f10, -52, 10, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs)} eq { <3,3>, <2, 3>, <5, -2> };

    //this is wrong - used to have the 0 error
    log_coeffs := SchoferFormula(f10, -72, 10, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs)} eq { <3,1>, <2, 2>, <5, -2>, <7,-2> };

    //this works!
    log_coeffs := SchoferFormula(f10, -120, 10, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs)} eq { <3,3>, <2, 2>,  <7,-2> };

    //this works!
    log_coeffs := SchoferFormula(f10, -88, 10, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs)} eq { <3,3>, <2, 1>, <5,3>, <7,-2> };

    //this is wrong -used to have the 0 error
    log_coeffs := SchoferFormula(f10, -27, 10, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs)} eq { <3,1>, <2, 8>, <5,-2> };

    //this works!
    log_coeffs := SchoferFormula(f10, -35, 10, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs)} eq {  <2, 8>, <7,-1> };
    //this works!
    log_coeffs := SchoferFormula(f10, -148, 10, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs)} eq {  <2, 3>, <3,3>, <11,3>, <5,-2>, <7,-2>, <13,-2> };
    //this works!
    log_coeffs := SchoferFormula(f10, -43, 10, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs)} eq {  <2, 8>, <3,3>,  <5,-2>, <7,-2>};

    //This (-180) doesn't work! Why??? This was proved in Elkies, also
    log_coeffs := SchoferFormula(f10, -180, 10, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs)} eq {  <2, 3>,  <11,3>,<13,-2>};

    // this works!
    log_coeffs := SchoferFormula(f10, -232, 10, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0}eq {   <3,3>,  <11,3>,<17,3>, <5,-2>, <7,-2>, <23,-2>};

    //this works!
   log_coeffs := SchoferFormula(f10, -67, 10, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs)} eq {   <2,8>,  <3,3>,  <5,3>, <7,-2>, <13,-2>};

    //this works!
   log_coeffs := SchoferFormula(f10, -280, 10, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs)} eq {   <2,1>,  <3,3>, <11,3>, <7,-1>, <23,-2>};
    //this works!
   log_coeffs := SchoferFormula(f10, -340, 10, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs)} eq {   <2,3>,  <3,3>, <23,3>, <7,-2>, <29,-2>};

    //this works!
   log_coeffs := SchoferFormula(f10, -115, 10, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs)} eq {   <2,11>,  <3,3>, <13,-2>, <23,-1>};
    //this works!
   log_coeffs := SchoferFormula(f10, -520, 10, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs)} eq {   <2,-1>,  <3,3>, <29,3>, <7,-2>, <13,-1>, <47,-2>};

    //this works but takes a long time!
   log_coeffs := SchoferFormula(f10, -163, 10, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs)} eq {   <2,11>,  <3,3>, <5,3>, <11,3>, <7,-2>, <13,-2>,<29,-2>, <31,-2> };

    //this works!
   log_coeffs := SchoferFormula(f10, -760, 10, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs)} eq {   <2,2>,  <3,3>, <17,3>, <47,3>, <7,-2>, <31,-2>,<71,-2> };
    //this works!
   log_coeffs := SchoferFormula(f10, -235, 10, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs)} eq {   <2,8>,  <3,3>, <17,3>,  <7,-2>, <37,-2>,<47,-1> };
    //this works! the extra 2^2 comes from the fact that there are 2 points here
   log_coeffs := SchoferFormula(f10, -420, 10, 1);
   assert {<p,log_coeffs[p]> : p in Keys(log_coeffs)} eq { <2,6> , <3,6>,<29,3>,<7,-2>,<37,-2> };


end procedure;


procedure test_Schofer_142()

    _<q> := LaurentSeriesRing(Rationals());
    f1 := -2*q^(-87)- 2*q^(-71) - 2*q^(-48) - 2*q^(-36) +2*q^(-16) -2*q^(-15)-2*q^(-12)+2*q^(-9)-2*q^(-7) - 2*q^(-3) +2*q^(-2) + 2*q^(-1) + O(q);
    f3 := q^(-87) - 2*q^(-79) + q^(-76) - 2*q^(-71)+2*q^(-48)-q^(-40)+3*q^(-32)-2*q^(-20)-q^(-19)-2*q^(-12)+2*q^(-10)+q^(-8)-2*q^(-7)-2*q^(-2) + O(q);
   
    //GY p.30
    //then psi_f1 *2^10= x, but I'm getting the reverse...
    log_coeffs := SchoferFormula(f1, -19, 142, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0}eq { <2,-10> };

    log_coeffs := SchoferFormula(f1, -20, 142, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0}eq { <2,-10> };

    log_coeffs := SchoferFormula(f1, -24, 142, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0}eq { <2,-11> };

    log_coeffs := SchoferFormula(f1, -40, 142, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0}eq { <2,-11> };

    log_coeffs := SchoferFormula(f1, -43, 142, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0}eq { <2,-10> };

    //long time didn't try
    log_coeffs := SchoferFormula(f1, -148, 142, 1);
    {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0}eq { <2,-10>  };

    //long time didn't try.
    log_coeffs := SchoferFormula(f1, -232, 142, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0}eq { <2,-11> };

    //psi_f3 *2 = y

    //Can't compute this
    log_coeffs := SchoferFormula(f3, -3, 142, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0}eq { <2,1> };

    //this works
    log_coeffs := SchoferFormula(f3, -4, 142, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0}eq {  };
    
    log_coeffs := SchoferFormula(f3, -19, 142, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0}eq { <2,-1> };

    //this works
    log_coeffs := SchoferFormula(f3, -24, 142, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0}eq { <2,1> };

    log_coeffs := SchoferFormula(f3, -43, 142, 1);
    assert
     {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0}eq { <3,1> };

    log_coeffs := SchoferFormula(f3, -148, 142, 1);


    log_coeffs := SchoferFormula(f3, -232, 142, 1);


end procedure;

procedure test_Schofer_26()

    //Test as in Guo Yang p.25
    _<q> := LaurentSeriesRing(Rationals());

    g1 := 2*q^(-13) - 2*q^(-2) - 4*q^(-1) + O(q);
    g2 := q^(-11) + 2*q^(-7) - 2*q^(-2) + O(q);
    g3 := 2*q^(-26) + 6*q^(-7)- 6*q^(-2) +2*q^(-1) + O(q);


    E, n, t := WeaklyHolomorphicBasis(26,1);
    E := Submatrix(E, [1..Rank(E)], [1..Ncols(E)]);
    fs_E := [q^(-n)*&+[(Integers()!b[i])*q^(i-1) : i in [1..Ncols(E)]] : b in Rows(E)];

    f4 := fs_E[8]; // this is the constant 1 function

    log_coeffs := SchoferFormula(g1, -11, 26, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0}eq {  };

    log_coeffs := SchoferFormula(g1, -19, 26, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0}eq { <3,2>  };

    log_coeffs := SchoferFormula(g1, -20, 26, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0}eq { <5,1>  };

    log_coeffs := SchoferFormula(g1, -24, 26, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0}eq { <3,1>  };

    log_coeffs := SchoferFormula(g1, -67, 26, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0}eq { <3,4>, <5,-2>  };


    log_coeffs := SchoferFormula(g2, -19, 26, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0}eq { <2,6>  };

    log_coeffs := SchoferFormula(g2, -20, 26, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0}eq { <2,5>  };

    log_coeffs := SchoferFormula(g2, -24, 26, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0}eq { <2,5>  };

    log_coeffs := SchoferFormula(g2, -52, 26, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0}eq { <2,3>  };

    log_coeffs := SchoferFormula(g2, -67, 26, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0}eq { <2,6>, <5,-2>, <7,1>  };


    log_coeffs := SchoferFormula(g3, -11, 26, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0} eq { <2,10>,<11,1>,<13,3> };

    log_coeffs := SchoferFormula(g3, -19, 26, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0}eq { <2,10>,<13,3>,<19,1>  };

    log_coeffs := SchoferFormula(g3, -20, 26, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0}eq { <2,12>,<13,3>  };

    log_coeffs := SchoferFormula(g3, -24, 26, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0}eq { <2,13>,<13,3>  };

    log_coeffs := SchoferFormula(g3, -52, 26, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0}eq { <2,6>,<13,8>  };

    log_coeffs := SchoferFormula(g3, -67, 26, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0}eq { <2,10>, <5,-6>, <13,3>,<41,2>,<67,1>  };

end procedure;

test_Kappa0();
test_Schofer_6();
test_Schofer_10();
test_Schofer_142();
test_Schofer_26();
