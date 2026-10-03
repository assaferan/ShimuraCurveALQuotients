// A cover X_0(D,N)/W is the fibre product, over the star curve, of its index-2 Atkin-Lehner
// double covers.  This test exercises that construction on committed models at X_0(6,23):
//
//   * the genus-5 TOP curve, which has no model because it is not subhyperelliptic (no genus-0
//     Atkin-Lehner quotient), comes out of three double covers and matches the Eichler-Selberg
//     trace formula on the 6-new space at four primes over two degrees;
//   * W={1,6}, which the pipeline DID build and stores as a pair, is reproduced by the same code
//     (a positive control: the construction must agree with the existing route where both apply);
//   * at X_0(21,2) the committed double covers come from runs with different Hauptmodul
//     normalisations, and their fibre product has the right genus but the WRONG point counts, so
//     the stage must refuse it (a negative control: the genus alone is not an acceptance test).
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
models_21_2 := eval (Read("data/models/models_21_2.m") cat "\nreturn models;");
poly_6_23 := hyperelliptic_polys(models_6_23);
poly_21_2 := hyperelliptic_polys(models_21_2);

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

procedure test_fibre_product_6_23(poly)
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

    // POSITIVE CONTROL: W={1,6} (genus 3), built by the pipeline as a pair, rebuilt here
    X6 := rep{Y : Y in curves | Y`D eq 6 and Y`N eq 23 and Y`W eq {Integers()|1,6}};
    ok6, gens6 := FibreProductGenerators({Integers()|1,6},
                     {U : U in AtkinLehnerDoubleCoversOver({Integers()|1,6}, 6, 23) | IsDefined(poly, U)}, 6, 23);
    assert ok6 and #gens6 eq 2;
    K6 := FibreProductFunctionField([poly[U] : U in gens6]);
    assert Genus(K6) eq X6`g and X6`g eq 3;

    // THE STAGE ITSELF: it must fill the absent top curve with a genus-5 curve given by 3 equations
    all_eqns, labels := pipeline_shaped(poly, 6, 23, curves);
    itop := rep{i : i in labels | curves[i]`W eq {Integers()|1}};
    assert not IsDefined(all_eqns, itop);
    all_ws := AssociativeArray();
    all_eqns, all_ws := EquationsByFibreProduct(all_eqns, all_ws, curves : NPrimes := 2);
    assert IsDefined(all_eqns, itop) and IsDefined(all_eqns[itop], 0);
    assert #DefiningPolynomials(all_eqns[itop][0]) eq 3;
end procedure;

procedure test_fibre_product_refuses_mixed_coordinates_21_2(poly)
    curves := GetHyperellipticCandidates();
    all_eqns, labels := pipeline_shaped(poly, 21, 2, curves);
    i14 := rep{i : i in labels | curves[i]`W eq {Integers()|1,14}};
    assert not IsDefined(all_eqns, i14);
    // the committed index-2 covers above W={1,14} do give a compositum of the right genus ...
    ups := AtkinLehnerDoubleCoversOver({Integers()|1,14}, 21, 2);
    ok, gens := FibreProductGenerators({Integers()|1,14}, {U : U in ups | IsDefined(poly, U)}, 21, 2);
    assert ok;
    assert Genus(FibreProductFunctionField([poly[U] : U in gens])) eq curves[i14]`g;
    // ... but it is not X_0(21,2)/w_14, and the stage must leave the key alone
    all_ws := AssociativeArray();
    all_eqns, all_ws := EquationsByFibreProduct(all_eqns, all_ws, curves : NPrimes := 2);
    assert not IsDefined(all_eqns, i14);
end procedure;

test_fibre_product_6_23(poly_6_23);
test_fibre_product_refuses_mixed_coordinates_21_2(poly_21_2);
