// tests/MolinaEquations.m
//
// Published equations of hyperelliptic Shimura curves and their Atkin-Lehner quotients, checked
// against this repository.  They are an outside source of a kind we did not have before: Molina
// derives them from Cerednik-Drinfeld uniformisation at a prime dividing D, together with the
// specialisation of Heegner points and Gross-Zagier's formula for their endomorphism rings, where
// every equation of ours comes from Borcherds products.  The two constructions share no step.
//
// SOURCE.  S. Molina, "Equations of hyperelliptic Shimura curves", arXiv:1004.3675v1 (2010),
// Section 10.4, Table 2, transcribed verbatim below.  He notes these equations had previously only
// been conjectured, by Kurihara ("On p-adic Poincare series and Shimura curves", Internat. J. Math.
// 5 (1994)).
//
// WHAT IS CHECKED.
//  [1] Each published curve has the genus the genus formula predicts, and its numbers of points
//      over F_p and F_{p^2} agree with the trace formula at every odd prime of good reduction
//      below 40.  Neither side knows about the other, so agreement tests both.
//  [2] Where we store a model of the same quotient, the two curves are isomorphic over Q.  The
//      coordinates differ -- ours is written over a hauptmodul of the star curve and his is not --
//      so the comparison is by isomorphism, not by coefficients.
//  [3] X_0^85(1)/<w_5> and X_0^85(1)/<w_85> are in his table and we hold no model for them, so
//      those two rows are checked against the trace formula only.  They are equations this
//      repository does not otherwise have.

printf "Testing the published equations of Molina...";

P<x> := PolynomialRing(Rationals());

// <D, W, published genus, y^2 = f(x)>, exactly as printed in Table 2
molina := [* <39, {Integers()|1},     3, -(7*x^4+79*x^3+311*x^2+497*x+277)*(x^4+9*x^3+29*x^2+39*x+19)>,
             <55, {Integers()|1},     3, -(3*x^4-32*x^3+130*x^2-237*x+163)*(x^4-8*x^3+34*x^2-83*x+81)>,
             <35, {Integers()|1,5},   2, -x*(9*x+4)*(4*x+1)*(172*x^3+176*x^2+60*x+7)>,
             <51, {Integers()|1,17},  2, -x*(7*x^3+52*x^2+116*x+68)*(x-1)*(x+3)>,
             <57, {Integers()|1,3},   2, -(x-9)*(x^3-19*x^2+119*x-249)*(7*x^2-104*x+388)>,
             <65, {Integers()|1,13},  2, -(x^2-3*x+1)*(7*x^4-3*x^3-32*x^2+25*x-5)>,
             <65, {Integers()|1,5},   2, -(x^2+7*x+9)*(7*x^4+81*x^3+319*x^2+508*x+268)>,
             <69, {Integers()|1,23},  2, -x*(x+4)*(4*x^4-16*x^3+11*x^2+10*x+3)>,
             <85, {Integers()|1,5},   2, -(3*x^2-41*x+133)*(x^4-23*x^3+183*x^2-556*x+412)>,
             <85, {Integers()|1,85},  2, (x^2-3*x+1)*(x^4+x^3-15*x^2+20*x-8)> *];

ncounts := 0;
nmodels := 0;
for e in molina do
    D := e[1]; W := e[2]; gpub := e[3]; f := e[4];
    C := HyperellipticCurve(f);

    // [1] genus, from the formula and from the published equation
    assert GenusShimuraCurveQuotient(D, 1, W) eq gpub;
    assert Genus(C) eq gpub;

    // [1] point counts against the trace formula
    X := CreateShimuraQuot(D, 1, W);
    X`g := gpub;
    for p in PrimesUpTo(40) do
        if p eq 2 or D mod p eq 0 then continue; end if;
        fp := ChangeRing(f, GF(p));
        if Discriminant(fp) eq 0 or Degree(fp) lt 2*gpub + 1 then continue; end if;
        for d in [1, 2] do
            assert #Points(ChangeRing(C, GF(p^d))) eq ComputePointsViaTrace(X, p, d);
            ncounts +:= 1;
        end for;
    end for;

    // [2] against our own model of the same quotient, where we have one
    fn := Sprintf("data/models/models_%o_1.m", D);
    if not FileExists(fn) then continue; end if;
    models := eval (Read(fn) cat "\nreturn models;");
    k := Sort(SetToSequence(W));
    if not IsDefined(models, k) or #models[k] eq 0 then continue; end if;
    ours := models[k][1];
    if Type(ours[2]) eq MonStgElt then continue; end if;       // stored as a fibre product
    g := ours[2] + (#ours ge 3 select ours[3]^2/4 else 0);
    assert IsIsomorphic(HyperellipticCurve(g), C);
    nmodels +:= 1;
end for;

// [3] the two rows at D = 85 are the ones we hold no model for
assert not FileExists("data/models/models_85_1.m");
assert #[e : e in molina | e[1] eq 85] eq 2;

assert ncounts eq 180 and nmodels eq 8;
printf " ok (%o published equations, %o point counts against the trace formula, %o isomorphic to ours)\n",
       #molina, ncounts, nmodels;
