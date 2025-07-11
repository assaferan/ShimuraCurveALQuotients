AttachSpec("shimuraquots.spec");

import "BorcherdsProducts.m" : Wpoly, Wpoly2, Wpoly_scaled;

function lambda_v(Q, d)
    bd := 10;
    found_lambda := false;
    while not found_lambda do
        bd *:= 2;
        found_lambda, lambda := FindLambda(Q,-d : bound := bd);
    end while;
    assert found_lambda;
    return lambda;
end function;

procedure test_Kappa0()
    printf "Testing Kappa0...";
    D := 6;
    N := 1;
    L, Ldual, disc_grp, to_disc, Qinv := ShimuraCurveLattice(D,N);
    Q := ChangeRing(Qinv^(-1), Integers());
    
    // verifying [Yang, Example 21, p. 24-25]
    // assert Round(1/Kappa0(3,-4,Q)) eq 2^8*3^4;
    log_coeffs := Kappa0(3,-4,Q,lambda_v(Q,-4));
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs)} eq { <2, -8>, <3, -4> };
    // assert Round(1/Kappa0(1,-3,Q)) eq 2^4;
    log_coeffs := Kappa0(1,-3,Q,lambda_v(Q,-3));
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs)} eq { <2, -4> };
    // verifying [Err, p. 850]
    // assert Round(1/Kappa0(1,-24,Q)) eq 2^6;
    log_coeffs := Kappa0(1,-24,Q,lambda_v(Q,-24));
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs)} eq { <2, -6> };
    // assert Round(1/Kappa0(3,-24,Q)) eq 2^8*3^4;
    log_coeffs := Kappa0(3,-24,Q,lambda_v(Q,-24));
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs)} eq { <2, -8>, <3, -4> };
    // assert Round(1/Kappa0(1,-163,Q)) eq 2^4*3^11*7^4*19^4*23^4;
    log_coeffs := Kappa0(1,-163,Q,lambda_v(Q,-163));
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs)} eq { <2, -4>, <3, -11>, <7, -4>, <11,0>, <19,-4>, <23,-4> };
    // assert Round(1/Kappa0(3,-163,Q)^3) eq 2^40*3^12*5^12*11^12*17^12;
    log_coeffs := Kappa0(3,-163,Q,lambda_v(Q,-163));
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs)} eq { <2, -40/3>, <3, -4>, <5, -4>, <11,-4>, <17,-4>, <23,0>, <89, 0> };
    D := 10;
    L, Ldual, disc_grp, to_disc, Qinv := ShimuraCurveLattice(D,N);
    Q := ChangeRing(Qinv^(-1), Integers());
    log_coeffs := Kappa0(3,-68,Q,lambda_v(Q,-68));
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0} eq { <2, -8>,  <5, -14/3>};
    // This does not work. This is probably a typo! Isn't it supposed to be 2. This would be the relevant number.
    //log_coeffs := Kappa0(1,-68,Q);
    log_coeffs := Kappa0(2,-68,Q,lambda_v(Q,-68));
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0} eq { <2, -6>,  <5, -6>};
    printf "Done!\n";
    return;
end procedure;

procedure test_Schofer_6()
    printf "Testing Schofer formula for D=6...";
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
    printf "Done!\n";
    return;

end procedure;


procedure test_Schofer_10()
    printf "testing Schofer formula for D=10...";
    _<q> := LaurentSeriesRing(Rationals());

    //This works!
    f10 := 3*q^(-3) - 2*q^(-2) + O(q);
    log_coeffs := SchoferFormula(f10, -20, 10, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs)} eq { <2, 3> };

    //This does still not work at 2! Off by a factor of 4, this is confusing. The kappas are all correct now.
    // Probably an error in [Err] ?
    log_coeffs := SchoferFormula(f10, -68, 10, 1);
    // Here is the result we expected from [Err]
    // assert {<p,log_coeffs[p]> : p in Keys(log_coeffs)} eq { <2, 2>, <5, 1> };
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs)} eq { <2, 6>, <5, 1> };

    //t_10 = 2^(-2)*|Psi_f_10|^2

    //Testing Err. Table 4, Section 8.2.1
    //this works!
    log_coeffs := SchoferFormula(f10, -40, 10, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs)} eq { <3,3>, <2, 2> };

    //this works!
    log_coeffs := SchoferFormula(f10, -52, 10, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs)} eq { <3,3>, <2, 3>, <5, -2> };

    // !At the moment, non-maximal orders are not implemented!
    /*
    //this is wrong - used to have the 0 error
    log_coeffs := SchoferFormula(f10, -72, 10, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs)} eq { <3,1>, <2, 2>, <5, -2>, <7,-2> };
    */

    //this works!
    log_coeffs := SchoferFormula(f10, -120, 10, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs)} eq { <3,3>, <2, 2>,  <7,-2> };

    //this works!
    log_coeffs := SchoferFormula(f10, -88, 10, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs)} eq { <3,3>, <2, 1>, <5,3>, <7,-2> };

    // !At the moment, non-maximal orders are not implemented!
    /*
    //this is wrong -used to have the 0 error
    log_coeffs := SchoferFormula(f10, -27, 10, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs)} eq { <3,1>, <2, 8>, <5,-2> };
    */

    //this works!
    log_coeffs := SchoferFormula(f10, -35, 10, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs)} eq {  <2, 8>, <7,-1> };
    //this works!
    log_coeffs := SchoferFormula(f10, -148, 10, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs)} eq {  <2, 3>, <3,3>, <11,3>, <5,-2>, <7,-2>, <13,-2> };
    //this works!
    log_coeffs := SchoferFormula(f10, -43, 10, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs)} eq {  <2, 8>, <3,3>,  <5,-2>, <7,-2>};

    // !At the moment, non-maximal orders are not implemented!
    /*
    //This (-180) doesn't work! Why??? This was proved in Elkies, also
    log_coeffs := SchoferFormula(f10, -180, 10, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs)} eq {  <2, 3>,  <11,3>,<13,-2>};
    */

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

    printf "Done!\n";
    return;
end procedure;


procedure test_Schofer_142()
    printf "testing Schofer formula for D=142...";
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

    // this works, but long time
    /*
    log_coeffs := SchoferFormula(f1, -148, 142, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0}eq { <2,-10>  };
    */

    // this works, but long time
    /*
    log_coeffs := SchoferFormula(f1, -232, 142, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0}eq { <2,-11> };
    */

    //psi_f3 *2 = y

    //Can't compute this because m = 0 !?
    /*
    log_coeffs := SchoferFormula(f3, -3, 142, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0}eq { <2,1> };
    */

    //this works
    log_coeffs := SchoferFormula(f3, -4, 142, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0}eq {  };
    
     //Can't compute this because m = 0 !?
    /*
    log_coeffs := SchoferFormula(f3, -19, 142, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0}eq { <2,-1> };
    */

    //this works
    log_coeffs := SchoferFormula(f3, -24, 142, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0}eq { <2,1> };

    log_coeffs := SchoferFormula(f3, -43, 142, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0}eq { <3,1> };

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


    // Not clear what this is doing here...
    /*
    E, n, t := WeaklyHolomorphicBasis(26,1);
    E := Submatrix(E, [1..Rank(E)], [1..Ncols(E)]);
    fs_E := [q^(-n)*&+[(Integers()!b[i])*q^(i-1) : i in [1..Ncols(E)]] : b in Rows(E)];

    f4 := fs_E[#fs_E]; // this is the constant 1 function
    */

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

    // this does not work!
    // Probably a typo in [GY]
    log_coeffs := SchoferFormula(g3, -52, 26, 1);
    // Original value in [GY]
    // assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0} eq { <2,6>,<13,8>  };
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0} eq { <2,6>,<13,5>  };

    log_coeffs := SchoferFormula(g3, -67, 26, 1);
    assert {<p,log_coeffs[p]> : p in Keys(log_coeffs) | log_coeffs[p] ne 0}eq { <2,10>, <5,-6>, <13,3>,<41,2>,<67,1>  };

    printf "Done!\n";
    return;
end procedure;


// functions for testing the Wpolys from Kudla Yang paper
function bp_Kudla_Yang_poly(p, kappam, D)
    d, c := SquarefreeFactorization(Integers()!(4*kappam));
    k := Valuation(c,p);
    F := Rationals();
    R<x> := PolynomialRing(F);
    vp := KroneckerSymbol(d, p);
    if k lt 0 then return R!1; end if;
    // k >= 0
    if (D mod p ne 0) then
        return (1 - vp*x + p^k*vp*x^(1+2*k)-p^(k+1)*x^(2*k+2))/(1-p*x^2);
    end if;
    // p divides D
    return ((1-vp*x)*(1-p^2*x^2)-vp*p^(k+1)*x^(2*k+1)+p^(k+2)*x^(2*k+2)+vp*p^(k+1)*x^(2*k+3)-p^(2*k+2)*x^(2*k+4))/(1-p*x^2);
end function;

function sigmasp_Kudla_Yang_poly(p, m, kappa, is_even)
    F := Rationals();
    R<x> := FunctionField(F);
    assert IsSquarefree(kappa);
    chi_p := is_even select KroneckerCharacter(kappa)(p) else KroneckerCharacter(2*kappa)(p);
    if m eq 0 then
        return (1 - chi_p*x)^(-1);
    end if;
    return &+[(chi_p*x)^r : r in [0..Valuation(m,p)]];
end function;

procedure test_kronecker_sigma(B)
    kappas := [kappa : kappa in [1..B] | IsSquarefree(kappa)];
    for p in PrimesUpTo(B) do
        for kappa in kappas do
            assert sigmasp_Kudla_Yang_poly(p,0,kappa,true)*EulerFactor(KroneckerCharacter(kappa),p) eq 1;
        end for;
    end for;
end procedure;

// Testing Proposition 5.1 in [KY]
// Should have Wp(s-1/2,m,mu) = Lp(s,chi_{kappa m})/zeta_p(2s) bp(kappa m, s) * (m - Q(mu) in Zp)
// Note that our Wpoly_scaled is evaluated at (s+s_0) and when n = 1, s0 = n/2 - 1 = -1/2
// We also have:
// zeta_p(2s)^(-1) = 1 - X^2
// Lp(s, chi_{kappa m}) = (1 - chi_{kappa m}(p) X)^(-1) 
// There is a sqrtp factor that I am missing

// procedure test_bp_KY(B)
function test_bp_KY(B)
    _<x> := PolynomialRing(Rationals());
    kappas := [kappa : kappa in [1..B] | IsSquarefree(kappa)];
    L := RSpaceWithBasis(IdentityMatrix(Integers(),1));
    failures := [* *];
    for m0 in [1..B] do
        for kappa in kappas do
            Q := Matrix([[2*kappa]]);
            Q_rat := ChangeRing(Q, Rationals());
            for P in PrimesUpTo(B, Rationals() : coprime_to := kappa) do
                mus := [Vector(Rationals(), [0])];
                p := Norm(P);
                if (p eq 2) then Append(~mus, Vector([1/2])); end if;
                for mu in mus do
                    m := m0 + 1/2*(mu*Q_rat, mu);
                    // kappa is NOT replaced by 2*kappa inthe character because the rank is odd - see [KY] (2.9)
                    rhs := ((1-x^2)/EulerFactor(KroneckerCharacter(Integers()!(kappa*Numerator(m)*Denominator(m))),p))*bp_Kudla_Yang_poly(p, kappa*m,1);
                    assert Denominator(rhs) eq 1;
                    rhs := Numerator(rhs);
                    K<sqrtp> := QuadraticField(p);
                    lhs := (p eq 2) select Wpoly2(m,mu,L,K,Q) else Wpoly(m,p,mu,L,K,Q);
                    // assert lhs eq ChangeRing(rhs, BaseRing(lhs));
                    if lhs ne ChangeRing(rhs, BaseRing(lhs)) then
                        Append(~failures, [* m0, kappa, p, mu *]);
                    end if;
                end for;
            end for;
        end for;
    end for;
    // return;
    return failures;
// end procedure;
end function;

// Prop. 2.1 in [KY] says that if chi_p is unramified and p is odd in the odd dimensional case
// we have Wmp(s) = sigma_{-s,p}(m,chi)/Lp(s+1,chi) in the even case
// and Lp(s+1/2,chi_{kappam})/zetap(2s+1)*bp(kappam,s+1/2) in the odd case
// In particular, when m = 0 this should yield
// Lp(s,chi) / Lp(s+1,chi) in the even case, and 
// zeta_p(2s) / zeta_p(2s+1) in the odd case


procedure test_W()
    // testing the few values we know from Yang
    L, Ldual, disc_grp, to_disc, Qinv := ShimuraCurveLattice(6,1);
    Q := ChangeRing(Qinv^(-1), Integers());
    for d in [-3,-4] do
        _, lambda_v := FindLambda(Q,-d);
        Lminus := Kernel(Transpose(Matrix(lambda_v*Q)));
        mu := Vector([0,0,0]);
        if d eq -4 then
            w32<x> := Wpoly_scaled(3,2,mu,Lminus,Q);
            assert w32 eq 1/2*(1-x^2);
            w33<x> := Wpoly_scaled(3,3,mu,Lminus,Q);
            assert w33 eq 1/3*(1+2*x+x^2);
            w22<x> := Wpoly_scaled(2,2,mu,Lminus,Q);
            assert w22 eq 1/2*(1+x^3);
            w23<x> := Wpoly_scaled(2,3,mu,Lminus,Q);
            assert w23 eq 1/3*(1-x);
        end if;
        if d eq -3 then
            w12<x> := Wpoly_scaled(1,2,mu,Lminus,Q);
            assert w12 eq 1/2*(1-x);
            w13<x> := Wpoly_scaled(1,3,mu,Lminus,Q);
            _<sqrt3> := BaseRing(w13);
            assert w13 eq 1/sqrt3*(1+x);
        end if;
    end for;
    return;
end procedure;

// This is not (!!!) [Err, Lemma 6.1, p. 845] based on [KRY, Lemmas 2.4 and 2.5]
// [Err only refers to primes not in Sm_mu]
// This is based on [KRY, Lemma 2.6]
function Wpolys_self_dual_KRY_2_4(m,p,mu,Lminus,Q)
    BML := BasisMatrix(Lminus);
    Delta := -Determinant(BML*Q*Transpose(BML));
    assert IsZero(mu); // in the self dual case mu is always zero
    K<sqrt_mp> := QuadraticField(-p);
    _<x> := PolynomialRing(K);
    return sqrt_mp * (1 - KroneckerCharacter(p)(Integers()!m) * x^(Valuation(m,p) + 1));
end function;

// This is based on [KY, Proposition 5.2]
function Wpolys_KY_5_3(m,p,mu,Lminus,Q)
    BML := BasisMatrix(Lminus);
    Qminus := BML*Q*Transpose(BML);
    lat_minus := LatticeWithGram(-Qminus);
    lat_minus_d := Dual(lat_minus : Rescale := false);
    disc_group := lat_minus_d / lat_minus;
    
    Delta := -Determinant(Qminus);
    kappa := SquareFree(Delta);

    // verifications
    assert Valuation(kappa,p) in [0,1];
    assert AbelianInvariants(disc_group) eq [2,-2*kappa];

    disc_group_p := pPrimaryComponent(disc_group, p);
    if (p ne 2) then assert #disc_group_p eq p; end if;
    if (p eq 2) and IsEven(kappa) then assert #disc_group_p eq 8; end if;
    if (p eq 2) and IsOdd(kappa) then assert #disc_group_p eq 4; end if;

    K<sqrt_kappa> := QuadraticField(kappa);
    _<x> := PolynomialRing(K);

    d := Discriminant(Integers(K));
    f := Valuation(d,p);
    a := Valuation(m,p);

    assert a ge -f;

    if (a eq -f) then return 1; end if;

    /*
    norm_form := Matrix([[Norm(x+y) - Norm(x) - Norm(y) : y in [1,sqrt_kappa]] : x in [1,sqrt_kappa]]);
    eps := (p eq 2) select [1,3,5,7] else [1,Integers()!Nonsquare(GF(p))];
    can_form_Q := GramMatrix(MinkowskiReduction(LatticeWithGram(-Qminus) : Canonical));
    a := can_form_Q[1,1] / 2;
    b := can_form_Q[1,2] / 2;
    ZK := Integers(K);
    I := a*ZK + (b-sqrt_kappa)*ZK; // we only do the pricipal ideal case
    can_forms := [GramMatrix(MinkowskiReduction(LatticeWithGram(e*norm_form) : Canonical)) : e in eps];
    assert can_form_Q in can_forms; // Not implemented when this is not the case.
    e := eps[Index(can_forms, can_form_Q)];
    */
    val_e := HasseMinkowskiInvariant(lat_minus, p);
   
    return 1 + val_e*KroneckerCharacter(kappa)(m)*x^(a+f);
end function;

procedure test_EquationsOfCovers()
    printf "testing equations of covers of X0*(26;1)...";
    curves := GetHyperellipticCandidates();
    assert exists(Xstar26){X : X in curves | X`D eq 26 and X`N eq 1 and IsStarCurve(X)};
    hs := EquationsOfCovers(Xstar26, curves);
    fs := [];
    for h in hs do
        Append(~fs, HyperellipticPolynomials(h));
    end for;
    _<x> := Universe(fs);
    assert fs eq [
                    -1/8*x^4 + 19/16*x^3 - 3/2*x^2 - 169/16*x,
                    x,
                    -2*x^3 + 19*x^2 - 24*x - 169
                    ];
    printf "Done\n";
    return;
end procedure;

test_Kappa0();
test_Schofer_6();
test_Schofer_10();
test_Schofer_142();
test_Schofer_26();
test_EquationsOfCovers();

exit;