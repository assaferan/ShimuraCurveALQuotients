// A cover X_0(D,N)/W is the fibre product, over the star curve, of its index-2 Atkin-Lehner
// double covers.  This test exercises that construction on committed models at X_0(6,23):
//
//   * the genus-5 TOP curve, which has no model because it is not subhyperelliptic (no genus-0
//     Atkin-Lehner quotient), comes out of three double covers and matches the Eichler-Selberg
//     trace formula on the 6-new space at four primes over two degrees;
//   * W={1,6}, which the pipeline DID build and stores as a pair of equations, is rebuilt by the
//     same code; the compositum matches the trace formula at four primes AND is compared with the
//     stored pair itself -- each stored equation y^2 = F(t) defines the same quadratic extension
//     of Q(t) as one of the three the rebuilt curve contains, so the two curves have the same
//     function field.  The same factors with one of them twisted by -1 give a curve of the same
//     genus that the trace formula rejects (the stage refuses a compositum the point counts
//     reject);
//   * at X_0(21,2) the double covers come from runs with different Hauptmodul normalisations, and
//     their fibre product has the right genus but the wrong point counts, so the stage must refuse
//     it (the genus alone does not accept a compositum).  Those polynomials are written out in
//     this file rather than read from the model file, so that the test keeps its meaning when that
//     file is regenerated;
//   * the triple check, which multiplies the three double covers above a common subgroup and asks
//     whether the product is a square: all seven covers of X_0(6,23) fall in one class, one cover
//     rescaled by 3 is caught, and six of the seven triples at X_0(21,2) fail.
//
// EXTERNAL SOURCE for the point counts: ComputePointsViaTrace, the codebase's own evaluation of
// the Eichler-Selberg trace formula on the W-fixed part of the D-new space, independent of every
// model-building stage.  Genera come from the Shimura-curve genus formula (GenusShimuraCurve).
//
// ⚠ The model files are evaluated here at TOP LEVEL, as every other test does.  Evaluating one
// inside a procedure body crashes Magma 2.29-7 when the test itself is run through run_tests.m,
// which evaluates the whole file as a string; the same file runs cleanly outside the harness, so
// the fault is invisible to a direct run.

_ := ClassNumberLU(-4);
P<x> := PolynomialRing(Rationals());

function hyperelliptic_polys(models)
    poly := AssociativeArray();
    for k in Keys(models) do
        if #models[k] gt 0 and Type(models[k][1][2]) ne MonStgElt then poly[Set(k)] := models[k][1][2]; end if;
    end for;
    return poly;
end function;

models_6_23 := eval (Read("data/models/models_6_23.m") cat "\nreturn models;");
models_14_19 := eval (Read("data/models/models_14_19.m") cat "\nreturn models;");
poly_6_23 := hyperelliptic_polys(models_6_23);

// The double covers of X_0(21,2)*, copied from data/models/models_21_2.m at commit faa2af2 (the
// first entry under each key) and written out here so that the test keeps its meaning if that file
// is ever regenerated in a single run.  The file mixes two Hauptmodul normalisations: the three
// covers above w_7 are in one coordinate, and the other four are each out of step with them, so of
// the seven triples of index-2 covers only {[1,2,7,14], [1,3,7,21], [1,6,7,42]} has a square
// product.  This is the mixed set the review of issue #99 reported and the one the stage refuses.
poly_21_2 := AssociativeArray();
poly_21_2[{Integers()|1,2}]          := P![ -4/9, -16/9, -23/9, -14/9, -7/9 ];
poly_21_2[{Integers()|1,3}]          := P![ 189/1024, 0, -1611/1024, 0, -801/1024, 0, -81/1024 ];
poly_21_2[{Integers()|1,6}]          := P![ 189/16, 0, 153/32, 0, -171/256, 0, -81/256 ];
poly_21_2[{Integers()|1,7}]          := P![ -3, 0, -1 ];
poly_21_2[{Integers()|1,14}]         := P![ 189/64, 0, 1053/128, 0, 5589/1024, 0, 189/1024 ];
poly_21_2[{Integers()|1,21}]         := P![ -7/256, 0, 31/128, 0, 9/256 ];
poly_21_2[{Integers()|1,42}]         := P![ -7/16, 0, -1/32, 0, 9/256 ];
poly_21_2[{Integers()|1,2,3,6}]      := P![ -4/9, 32/9, -71/9, 28/9 ];
poly_21_2[{Integers()|1,2,7,14}]     := P![ 0, 81/4, -81/4 ];
poly_21_2[{Integers()|1,2,21,42}]    := P![ -7/16, 3/32, 81/256 ];
poly_21_2[{Integers()|1,3,7,21}]     := P![ -3, 3 ];
poly_21_2[{Integers()|1,3,14,42}]    := P![ 0, 7/4, -1/8, -9/64 ];
poly_21_2[{Integers()|1,6,7,42}]     := P![ 0, -3 ];
poly_21_2[{Integers()|1,6,14,21}]    := P![ 0, 7/64, 31/32, -9/64 ];

// the two equations the pipeline stored for W={1,6}, as double covers y^2 = F(t) of the star line
stored_6_23 := CRVDoubleCovers(models_6_23[[Integers()|1,6]][1][3]);

// y^2 = a and y^2 = b are the same extension of Q(t) exactly when a*b is a square there: every
// irreducible factor to an even power, and the constant left over a square in Q
function same_quadratic_subfield(a, b)
    if a eq 0 or b eq 0 then return false; end if;
    fa, u := Factorization(a * b);
    if exists{v : v in fa | IsOdd(v[2])} then return false; end if;
    return IsSquare(u);
end function;

// every hyperelliptic committed model of one base over a single base label, pipeline-shaped
function pipeline_shaped(poly, D, N, curves)
    all_eqns := AssociativeArray();
    labels := [i : i in [1..#curves] | curves[i]`D eq D and curves[i]`N eq N];
    for i in labels do
        if IsDefined(poly, curves[i]`W) then
            all_eqns[i] := AssociativeArray();
            all_eqns[i][0] := HyperellipticCurve(poly[curves[i]`W]);
        end if;
    end for;
    return all_eqns, labels;
end function;

procedure test_degenerate_factor_sets()
    printf "  degenerate factor sets (a subset product that is a constant times a square)...";
    P<x> := PolynomialRing(Rationals());
    f1 := x^3 - x + 1; f2 := x^2 + 3;
    assert IsDegenerateFactorSet([f1, -3*f1]);            // product -3 f1^2: SquarefreePart would say f1
    assert IsDegenerateFactorSet([f1, f2, f1*f2]);        // product of all three is (f1 f2)^2
    assert IsDegenerateFactorSet([f1, 4*f1]);             // product 4 f1^2, a square: Magma would stop
    assert not IsDegenerateFactorSet([f1, f2]);
    assert not IsDegenerateFactorSet([f1, f2, f1*(x - 5)]);
    printf " ok\n";
end procedure;

procedure test_h_term_completes_the_square()
    // A stored double cover y^2 + h y = f with h != 0 enters the fibre product as y'^2 = f + h^2/4.
    // X_0^14(19), W = {1,2,19,38}: the one stored equation of the base with a nonzero h.  The function
    // field of y^2 + h y - f and that of y'^2 - (f + h^2/4) must have the same genus and the same
    // point counts over F_5 and F_7 (they are the same field).
    printf "  the h term: y^2 + h y = f enters as y^2 = f + h^2/4 (X0^14(19), W={1,2,19,38})...";
    e := models_14_19[[Integers()|1,2,19,38]][1];
    assert #e eq 3 and e[3] ne 0;
    f := e[2]; h := e[3];
    assert Genus(HyperellipticCurve(f, h)) eq e[1] and Genus(HyperellipticCurve(f + h^2/4)) eq e[1];
    checked := 0;
    for p in [5, 7, 11, 13, 17] do
        Kp := RationalFunctionField(GF(p)); R<Y> := PolynomialRing(Kp);
        L1 := FunctionField(Y^2 + Evaluate(PolynomialRing(GF(p))!h, Kp.1)*Y - Evaluate(PolynomialRing(GF(p))!f, Kp.1));
        L2 := FunctionField(Y^2 - Evaluate(PolynomialRing(GF(p))!(f + h^2/4), Kp.1));
        if Genus(L1) ne e[1] then continue; end if;         // bad reduction of this model at p
        assert Genus(L2) eq Genus(L1);
        assert [#Places(L1, d) : d in [1..2]] eq [#Places(L2, d) : d in [1..2]];
        checked +:= 1;
    end for;
    assert checked ge 2;
    printf " ok (%o good primes)\n", checked;
end procedure;

procedure test_fibre_product_6_23(poly, stored)
    curves := GetHyperellipticCandidates();

    // the seven index-2 subgroups above W={1}, all with hyperelliptic models
    ups := AtkinLehnerDoubleCoversOver({Integers()|1}, 6, 23);
    assert #ups eq 7;
    avail := {U : U in ups | IsDefined(poly, U)};
    assert #avail eq 7;

    // THE TOP CURVE, genus 5
    Xtop := rep{Y : Y in curves | Y`D eq 6 and Y`N eq 23 and Y`W eq {Integers()|1}};
    assert Xtop`g eq 5;
    ok, gens := FibreProductGenerators({Integers()|1}, avail, 6, 23);
    assert ok and #gens eq 3;
    fs := [poly[U] : U in gens];
    K := FibreProductFunctionField(fs);
    assert Genus(K) eq 5;
    // the Jacobian of a (Z/2)^3-cover of P^1 splits as the product over its seven double covers
    assert &+[Degree(poly[U]) le 2 select 0 else (Degree(poly[U]) - 1) div 2 : U in ups] eq 5;
    n_checked := 0;
    for p in [5, 7, 11, 13] do
        Kp := RationalFunctionField(GF(p)); L := Kp;
        for f in fs do
            R<Y> := PolynomialRing(L);
            L := FunctionField(Y^2 - L!Evaluate(PolynomialRing(GF(p))!f, Kp.1));
        end for;
        assert Genus(L) eq 5;
        cnt := [&+[e * #Places(L, e) : e in Divisors(d)] : d in [1..2]];
        exp := [ComputePointsViaTrace(Xtop, p, d) : d in [1..2]];
        assert cnt eq exp;
        n_checked +:= 1;
    end for;
    assert n_checked eq 4;
    C, eqns := FibreProductCurve(fs);
    assert #eqns eq 3;

    // W={1,6} (genus 3), built by the pipeline as a pair, is rebuilt here and checked against the
    // trace formula at four primes; then the same two factors with the first twisted by -1: same
    // genus, and the trace formula must reject it at some prime
    X6 := rep{Y : Y in curves | Y`D eq 6 and Y`N eq 23 and Y`W eq {Integers()|1,6}};
    ok6, gens6 := FibreProductGenerators({Integers()|1,6},
                     {U : U in AtkinLehnerDoubleCoversOver({Integers()|1,6}, 6, 23) | IsDefined(poly, U)}, 6, 23);
    assert ok6 and #gens6 eq 2;
    fs6 := [poly[U] : U in gens6];
    K6 := FibreProductFunctionField(fs6);
    assert Genus(K6) eq X6`g and X6`g eq 3;
    function counts_agree(fs, X, ps)
        agree := true; n := 0;
        for p in ps do
            Kp := RationalFunctionField(GF(p)); L := Kp;
            for f in fs do
                R<Y> := PolynomialRing(L);
                L := FunctionField(Y^2 - L!Evaluate(PolynomialRing(GF(p))!f, Kp.1));
            end for;
            if Genus(L) ne X`g then continue; end if;
            cnt := [&+[e * #Places(L, e) : e in Divisors(d)] : d in [1..2]];
            if cnt ne [ComputePointsViaTrace(X, p, d) : d in [1..2]] then agree := false; end if;
            n +:= 1;
        end for;
        return agree, n;
    end function;
    ok_pos, n_pos := counts_agree(fs6, X6, [5, 7, 11, 13]);
    assert ok_pos and n_pos eq 4;
    // and against the two equations the pipeline stored for this cover: the rebuilt curve has the
    // three quadratic subfields given by fs6[1], fs6[2] and their product, and each stored
    // equation must be one of them (two different ones, so the two compositums agree)
    // (control, 2026-10-13: scaling a stored equation by 4 still matches, as it must -- it is the
    // same extension -- while scaling it by 2 or by -1 matches nothing)
    assert #stored eq 2;
    subf := [fs6[1], fs6[2], fs6[1]*fs6[2]];
    matched := [];
    for f in stored do
        m := [j : j in [1..3] | same_quadratic_subfield(f, subf[j])];
        assert #m eq 1;
        Append(~matched, m[1]);
    end for;
    assert #Set(matched) eq 2;
    twisted := [-fs6[1], fs6[2]];
    assert Genus(FibreProductFunctionField(twisted)) eq 3;
    ok_neg, n_neg := counts_agree(twisted, X6, [5, 7, 11, 13]);
    assert (not ok_neg) and n_neg ge 2;

    // THE STAGE ITSELF: it must fill the absent top curve with a genus-5 curve given by 3 equations
    all_eqns, labels := pipeline_shaped(poly, 6, 23, curves);
    itop := rep{i : i in labels | curves[i]`W eq {Integers()|1}};
    assert not IsDefined(all_eqns, itop);
    all_ws := AssociativeArray();
    all_eqns, all_ws := EquationsByFibreProduct(all_eqns, all_ws, curves : NPrimes := 3);
    assert IsDefined(all_eqns, itop) and IsDefined(all_eqns[itop], 0);
    Ctop := all_eqns[itop][0];
    eqs := DefiningPolynomials(Ctop);
    assert #eqs eq 3;
    // the curve the stage wrote, read back from its equations alone (z = 1), against the trace
    // formula at three primes over two degrees
    A := Ambient(Ctop); Pt<t> := PolynomialRing(Rationals());
    fs_top := [];
    for i in [1..3] do
        F := eqs[i];
        // the i-th equation is y_i^2 - F_i(s, z): recover F_i(t, 1)
        yi := A.(2+i);
        Fi := yi^2 - F;
        assert Evaluate(Fi, [A.1, A.2] cat [0 : j in [1..3]]) eq Fi;
        Append(~fs_top, Evaluate(Fi, [t, 1] cat [0 : j in [1..3]]));
    end for;
    ok_top, n_top := counts_agree(fs_top, Xtop, [5, 7, 11]);
    assert ok_top and n_top eq 3;
    // and a rebased or back-filled key is never used as a factor
    all_eqns2, labels2 := pipeline_shaped(poly, 6, 23, curves);
    skip := {i : i in labels2 | IsDefined(all_eqns2, i)};
    all_eqns2, _ := EquationsByFibreProduct(all_eqns2, all_ws, curves : NPrimes := 2, Skip := skip);
    assert not IsDefined(all_eqns2, itop);
end procedure;

// The index-2 covers of one base, as the triple check wants them.
function index_two_covers(poly, D, N)
    DN := D * N;
    full := {Integers()| d : d in Divisors(DN) | GCD(d, DN div d) eq 1};
    covers := AssociativeArray();
    for U in Keys(poly) do
        if #U * 2 eq #full then covers[U] := poly[U]; end if;
    end for;
    return covers, full;
end function;

procedure test_triple_check_6_23(poly)
    // The seven double covers of X_0(6,23)* come from one pipeline run, so every triple of them
    // multiplies to a rational constant times a square and the check puts all seven in one class.
    printf "  the triple check on the seven double covers of X_0^6(23)...";
    covers := index_two_covers(poly, 6, 23);
    assert #Keys(covers) eq 7;
    failing, classes := AtkinLehnerTripleCheck(covers, 6, 23);
    assert #failing eq 0;
    assert #classes eq 1 and #classes[1] eq 7;
    printf " ok (no failing triple)\n";

    // Twisting one cover by 3 keeps its genus and its degree, and leaves the other six agreeing
    // with one another; the three triples through the twisted cover must fail.
    printf "  one cover of X_0^6(23) rescaled by 3 is detected...";
    R := Representative(Keys(covers));
    tw := AssociativeArray();
    for U in Keys(covers) do tw[U] := (U eq R select 3 else 1) * covers[U]; end for;
    failing, classes := AtkinLehnerTripleCheck(tw, 6, 23);
    assert #failing eq 3;
    assert #classes eq 2;
    assert exists{c : c in classes | c eq {R}} and exists{c : c in classes | #c eq 6};
    // a failing triple through the twisted cover: the other two share a class, so the twisted one
    // is the odd one and must not be multiplied with either
    tri := failing[1];
    mem := [tri[1], tri[2], tri[3]];
    assert R in mem;
    others := [U : U in mem | U ne R];
    assert DefinitelyDisagreeingFactors([R, others[1]], failing, classes);
    assert DefinitelyDisagreeingFactors([R, others[2]], failing, classes);
    // and the two untwisted members of the same triple are not accused of disagreeing
    assert not DefinitelyDisagreeingFactors(others, failing, classes);
    printf " ok\n";
end procedure;

procedure test_triple_check_mixed_coordinates_21_2(poly)
    printf "  the triple check on the mixed double covers of X_0^21(2)...";
    covers, full := index_two_covers(poly, 21, 2);
    assert #Keys(covers) eq 7;
    failing, classes := AtkinLehnerTripleCheck(covers, 21, 2);
    // one triple of the seven passes: the three covers above w_7
    assert #failing eq 6;
    w7 := {U : U in Keys(covers) | 7 in U};
    assert #w7 eq 3;
    assert exists{c : c in classes | c eq w7};
    assert #classes eq 5 and #[c : c in classes | #c eq 1] eq 4;

    // the three covers above W = {1, 14}, from any two of which the stage would build the cover,
    // are one of the failing triples, so they are not all in one coordinate
    ups := AtkinLehnerDoubleCoversOver({Integers()|1,14}, 21, 2);
    assert #ups eq 3;
    assert exists{tri : tri in failing | {tri[1], tri[2], tri[3]} eq ups};
    ok, gens := FibreProductGenerators({Integers()|1,14}, ups, 21, 2);
    assert ok and #gens eq 2;
    assert "inconclusive" in TripleCheckVerdict(gens, failing, classes);
    // The check does not say WHICH of the three is out of step: that would need the third cover to
    // share a class with one of the other two, and at this base the only class with more than one
    // member is the one above w_7, which no failing triple meets twice.  So the stage refuses this
    // cover on its point counts, not on the triple check.
    assert not DefinitelyDisagreeingFactors(gens, failing, classes);
    printf " ok (6 of 7 triples fail)\n";
end procedure;

procedure test_fibre_product_refuses_mixed_coordinates_21_2(poly)
    curves := GetHyperellipticCandidates();
    all_eqns, labels := pipeline_shaped(poly, 21, 2, curves);
    i14 := rep{i : i in labels | curves[i]`W eq {Integers()|1,14}};
    // The model file stores the genus-2 curve X_0^21(2)/w_14 (checked against Guo-Yang's equations
    // in tests/X0_21_2.m).  This test is about the stage, so that equation is taken out of the
    // input: the stage must not rebuild the cover from the committed double covers.
    if IsDefined(all_eqns, i14) then Remove(~all_eqns, i14); end if;
    // the committed index-2 covers above W={1,14} do give a compositum of the right genus ...
    ups := AtkinLehnerDoubleCoversOver({Integers()|1,14}, 21, 2);
    ok, gens := FibreProductGenerators({Integers()|1,14}, {U : U in ups | IsDefined(poly, U)}, 21, 2);
    assert ok;
    assert Genus(FibreProductFunctionField([poly[U] : U in gens])) eq curves[i14]`g;
    // ... but it is not X_0(21,2)/w_14 (its factors are in mixed Hauptmodul coordinates): the point
    // counts disagree with the trace formula at p = 11 (and agree at 13 and 23, which is why the
    // check runs over several good primes), and the stage must leave the key alone
    all_ws := AssociativeArray();
    all_eqns, all_ws := EquationsByFibreProduct(all_eqns, all_ws, curves : NPrimes := 2);
    assert not IsDefined(all_eqns, i14);
end procedure;

test_fibre_product_6_23(poly_6_23, stored_6_23);
test_fibre_product_refuses_mixed_coordinates_21_2(poly_21_2);

test_triple_check_6_23(poly_6_23);
test_triple_check_mixed_coordinates_21_2(poly_21_2);
test_degenerate_factor_sets();
test_h_term_completes_the_square();
