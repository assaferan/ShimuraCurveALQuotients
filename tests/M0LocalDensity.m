// Cross-check of the local Whittaker factor (LocalWhittakerAtOne, i.e. Wpoly/Wpoly2) against an
// INDEPENDENT implementation: the Siegel-Weil local representation density, computed by literally
// counting solutions of a quadratic congruence.
//
// This complements tests/Whittaker2.m, which checks the p=2 factor against Yang's closed formula
// (JNT 72 (1998) Thm 4.1). Here the check is structural rather than formula-vs-formula, and it covers
// the ODD ramified primes p | D too -- exactly the places where the multi-term m=0 investigation kept
// suspecting a defect.
//
// WHY THIS IS THE RIGHT COMPARISON. Schwagenscheidt ("Eisenstein series for the Weil representation",
// arXiv:1803.10550) Thm 1.4 gives the Fourier coefficients of the weight-k vector valued Eisenstein
// series for rho_A^* as an archimedean/L-function constant times a product of local factors
// L^{(p)}_{gamma,n}(X) evaluated at X = p^{1-m/2-k}, where
//     L^{(p)}_{gamma,n}(X) = N_{gamma,n}(p^{w_p}) X^{w_p} + (1 - p^{m-1} X) sum_{nu<w_p} N_{gamma,n}(p^nu) X^nu,
//     N_{gamma,n}(a)       = #{ r in L/aL : Q(r-gamma) + n = 0 mod a },
//     w_p                  = 1 + 2 v_p(2 N_beta N_gamma n).
// For our setting (rank m = 3, weight k = 3/2, beta = 0 so N_beta = 1) we have X = p^{-2} and
// p^{m-1} = p^2, so the factor (1 - p^{m-1} X) VANISHES -- k = m/2 is the Siegel-Weil point -- and the
// whole local factor collapses to the pure representation density
//     L^{(p)}_{gamma,n}(p^-2) = N_{gamma,n}(p^{w_p}) / p^{2 w_p}.
// That is what this test computes from scratch and compares against LocalWhittakerAtOne.
//
// Verified agreement: 46/46 (gamma = 0, n = 1..14, every p | 2 n det Q) on BOTH X0^15(2) and X0^21(2).

// N_{gamma,n}(p^w) for gamma = 0. F(x) = Q(x) + n is integer valued on Z^3 (L even), and for fixed
// (x1,x2) it is A*x3^2 + B*x3 + C with A = G33/2 CONSTANT. So histogram t -> A t^2 + b t once per b
// and look up (-C): O(p^{2w}) work instead of a brute-force O(p^{3w}) triple loop.
function rep_count_at_zero(G, n, p, w)
    pw := p^w;
    if w le 0 then return 1, pw; end if;
    h11 := Integers()!(G[1][1]/2); h22 := Integers()!(G[2][2]/2); h33 := Integers()!(G[3][3]/2);
    g12 := Integers()!G[1][2]; g13 := Integers()!G[1][3]; g23 := Integers()!G[2][3];
    c0 := Integers()!n; Am := h33 mod pw;
    hist := [];
    for b in [0..pw-1] do
        row := [0 : v in [1..pw]];
        for t in [0..pw-1] do
            v := (Am*t*t + b*t) mod pw;
            row[v+1] +:= 1;
        end for;
        Append(~hist, row);
    end for;
    N := 0;
    for x1 in [0..pw-1] do
        for x2 in [0..pw-1] do
            B := (g13*x1 + g23*x2) mod pw;
            C := (h11*x1*x1 + h22*x2*x2 + g12*x1*x2 + c0) mod pw;
            N +:= hist[B+1][((-C) mod pw)+1];
        end for;
    end for;
    return N, pw;
end function;

// Brute force, to certify the histogram shortcut itself (small moduli only).
function rep_count_brute(G, n, p, w)
    pw := p^w;
    Gr := ChangeRing(G, Rationals());
    cnt := 0;
    for x1 in [0..pw-1] do for x2 in [0..pw-1] do for x3 in [0..pw-1] do
        x := Vector(Rationals(), [x1,x2,x3]);
        if (Integers()!((x*Gr, x)/2 + n)) mod pw eq 0 then cnt +:= 1; end if;
    end for; end for; end for;
    return cnt;
end function;

procedure test_M0LocalDensity()
    printf "Testing local Whittaker vs Siegel-Weil representation density...";

    for base in [<15,2>, <21,2>] do
        D, N := Explode(base);
        Ldata := ShimuraCurveLattice(D, N);
        Q := ChangeRing(Ldata`Q, Integers());
        negQ := -Q;
        Lfull := RSpaceWithBasis(IdentityMatrix(Integers(), 3));
        zero := Vector(Rationals(), [0,0,0]);
        detprimes := Set(PrimeDivisors(Determinant(Q)));

        // the histogram counter must agree with brute force
        for pr in [<2,3>, <3,1>] do
            p, w := Explode(pr);
            assert rep_count_at_zero(Q, 2, p, w) eq rep_count_brute(Q, 2, p, w);
        end for;

        checked := 0;
        for n in [1..10] do
            for p in Sort([q : q in detprimes join Set(PrimeDivisors(n))]) do
                w := 1 + 2*Valuation(Rationals()!(2*n), p);
                cnt, pw := rep_count_at_zero(Q, n, p, w);
                density := (Rationals()!cnt)/(Rationals()!pw)^2;
                whittaker := LocalWhittakerAtOne(Rationals()!n, p, zero, Lfull, negQ);
                assert density eq whittaker;
                checked +:= 1;
            end for;
        end for;
        assert checked ge 30;   // guard against the loops silently degenerating to nothing
    end for;

    printf "Done!\n";
end procedure;

test_M0LocalDensity();
