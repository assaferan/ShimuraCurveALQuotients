// X_0^3390(1)/<w_5, w_6, w_113> IS HYPERELLIPTIC, PART ONE: THE COVER, ITS BRANCH LOCUS, AND THE TORSION OF THE STAR.
// About half a minute.
//
// Y = X_0^3390(1)/H, H = <w_5, w_6, w_113>, is the only genus-four quotient of X = X_0^3390(1) by an index-two subgroup
// of its Atkin-Lehner group W (order 16), and it is a double cover of the genus-two star X* = X/W.  By Riemann-Hurwitz
// the cover has two branch points; a point of the star is a branch point exactly when its orbit upstairs consists of
// fixed points of some w_m with m outside H, and Ogg's count of those fixed points by CM discriminant finds them only
// for w_3 (discriminant -3) and w_10 (discriminant -40), eight each: the branch points are the rational CM points
// P_-3 and P_-40 of the star.  Y is hyperelliptic if and only if P_-3 + P_-40 is a fibre of the hyperelliptic map of
// X* (paper/star-2002-proof.tex, Lemma 3), i.e. iff P_-3 + P_-40 - F is principal for a fibre F.  The rational torsion
// of J(X*) is trivial: the orders #J(X*)(F_p) = L_p(1), from the point counts over F_p and F_{p^2}, have gcd 1 over the
// good primes below 120.  So principal is the same as torsion, which tests/_offline/Hyperelliptic3390_1.m decides
// (paper/hyperelliptic-3390-1.tex).  This test also records the traces of Frobenius of the star used there to
// recognise the space W: T(49) has characteristic polynomial x^2 + 7x + 11 on it.

procedure test_hyperelliptic_3390_1_branch()
    printf "X_0^3390(1)/<w_5, w_6, w_113>: the unique genus-four quotient, branched over the star at P_-3 and P_-40, and J(X*)(Q) is torsion-free...";
    D := 3390; N := 1;
    W := { m : m in Divisors(D*N) | GCD(m, (D*N) div m) eq 1 };
    H := {1, 5, 6, 30, 113, 565, 678, 3390};
    assert #W eq 16 and H subset W and #H eq 8;
    assert GenusShimuraCurveQuotient(D, N, W) eq 2 and GenusShimuraCurveQuotient(D, N, H) eq 4;
    // the index-two subgroups are the kernels of the characters of W; H is the only one of genus four
    gens := [2, 3, 5, 113];
    for chi in [ [c1, c2, c3, c4] : c1, c2, c3, c4 in [0, 1] | [c1, c2, c3, c4] ne [0, 0, 0, 0] ] do
        H1 := { m : m in W | IsEven(&+([0] cat [ chi[i] : i in [1..4] | gens[i] in PrimeDivisors(m) ])) };
        assert (GenusShimuraCurveQuotient(D, N, H1) eq 4) eq (H1 eq H);
    end for;
    // the branch locus
    fixed := AssociativeArray();
    for m in W diff H do
        e := NumFixedPointsByCMOrder(D, N, m);
        for d in Keys(e) do if e[d] ne 0 then fixed[d] := <m, e[d]>; end if; end for;
        assert NumFixedPoints(D, N, m) eq &+([0] cat [ e[d] : d in Keys(e) ]);
    end for;
    assert Keys(fixed) eq {-3, -40} and fixed[-3] eq <3, 8> and fixed[-40] eq <10, 8>;
    Xstar := CreateShimuraQuot(D, N, W); Xstar`g := GenusShimuraCurveQuotient(D, N, W); Xstar`CurveID := 0;
    curves := GetQuotientsAndGenera([Xstar]);
    _ := exists(star){c : c in curves | IsStarCurve(c)};
    assert forall{d : d in [-3, -40, -43, -67, -123, -168] | DegreeOfFieldOfDefinitionOfCMPoint(star, d) eq 1};
    // the torsion of J(X*)(Q), and the Hecke data of the star at 7
    g := 0;
    for p in [q : q in PrimesInInterval(7, 120) | D mod q ne 0] do
        s1 := ComputePointsViaTrace(Xstar, p, 1); s2 := ComputePointsViaTrace(Xstar, p, 2);
        c1 := s1 - p - 1; c2 := (s2 - p^2 - 1 + c1^2) div 2;
        g := GCD(g, 1 + c1 + c2 + p*c1 + p^2);
        if p eq 7 then assert c1 eq 7 and c2 - 2*p eq 11; end if;
    end for;
    assert g eq 1;
    printf " Done!\n";
end procedure;

test_hyperelliptic_3390_1_branch();
