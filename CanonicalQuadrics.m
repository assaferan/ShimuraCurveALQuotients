// CanonicalQuadrics.m
//
// A pipeline stage that decides hyperellipticity from the canonical embedding, by Max Noether's
// count: the pairwise products of a basis of holomorphic differentials of a curve of genus g >= 3
// satisfy (g-1)(g-2)/2 linear relations when the curve is hyperelliptic (its canonical image is a
// rational normal curve of degree g-1) and (g-2)(g-3)/2 otherwise.  The count needs the
// differentials and their products, which the point-count filters never see, and it is exact.
//
// Two sources of differentials:
//   * D = 1: the holomorphic differentials of X_0(N)/W are the weight-2 cusp forms on Gamma_0(N)
//     fixed by the Atkin-Lehner involutions in W, as q-expansions; products of two are weight-4
//     forms, so coefficients past the Sturm bound of weight 4 determine every relation.
//   * D > 1: a model of the curve, from data/models/models_D_N.m: a stored y^2 + h y = f is
//     hyperelliptic outright; a stored list of equations (a fibre product of double covers of the
//     star line) or the fibre product of the stored double covers above W gives a function field
//     and the count runs there (IsHyperellipticFunctionField).  A fibre product is used only when
//     its factors pass the coordinate check (AtkinLehnerTripleCheck), it has the genus the formula
//     predicts and its point counts agree with the trace formula over F_p and F_{p^2} at three good
//     primes.  At N > 1 the stored double covers rest on the m = 0 term of Schofer's formula, and the
//     verdict's label says so.

// The W-fixed weight-2 cusp forms on Gamma_0(N) as q-expansions to precision prec.
function fixed_cusp_forms(N, W, prec)
    S := CuspForms(Gamma0(N), 2);
    if Dimension(S) eq 0 then return [], 0; end if;
    V := VectorSpace(Rationals(), Dimension(S));
    fixed := V;
    for q in W do
        if q eq 1 then continue; end if;
        fixed := fixed meet Kernel(AtkinLehnerOperator(S, q) - 1);
    end for;
    B := Basis(S);
    qexp := [PowerSeries(b, prec) : b in B];
    return [&+[v[i] * qexp[i] : i in [1..#B]] : v in Basis(fixed)], Dimension(fixed);
end function;

intrinsic HyperellipticByQExpansions(N::RngIntElt, W::SetEnum) -> BoolElt, RngIntElt, RngIntElt
{For X_0(N)/W of genus g >= 3: whether it is hyperelliptic, by Noether's count on the W-fixed
 weight-2 cusp forms; also returns the number of quadratic relations found and the dimension of
 the fixed space (which must be g).}
    idx := N * &*[1 + 1/p : p in PrimeDivisors(N)];
    prec := Ceiling(4 * idx / 12) + 10;                 // past the Sturm bound of weight 4
    fs, g := fixed_cusp_forms(N, W, prec);
    require g ge 3 : "the fixed space has dimension less than 3";
    prods := [fs[i] * fs[j] : j in [i..g], i in [1..g]];
    M := Matrix(Rationals(), [[Coefficient(f, n) : n in [0..prec-1]] : f in prods]);
    nrel := #prods - Rank(M);
    hyp := nrel eq (g-1)*(g-2) div 2;
    require hyp or nrel eq (g-2)*(g-3) div 2 : "the number of quadratic relations matches neither case";
    return hyp, nrel, g;
end intrinsic;

// The stored double covers of a base, index-2 subgroup -> f with h^2/4 added, and the entry under
// a key, from the model file; false when there is no file.
function stored_models(D, N, ModelDir)
    fn := Sprintf("%o/models_%o_%o.m", ModelDir, D, N);
    if not FileExists(fn) then return false, _, _; end if;
    P<x> := PolynomialRing(Rationals());
    models := eval (Read(fn) cat "\nreturn models;");
    DN := D*N;
    full := {Integers()| d : d in Divisors(DN) | GCD(d, DN div d) eq 1};
    covers := AssociativeArray();
    for k in Keys(models) do
        if #models[k] gt 0 and Type(models[k][1][2]) ne MonStgElt and #Set(k)*2 eq #full then
            e := models[k][1];
            covers[Set(k)] := e[2] + (#e ge 3 select e[3]^2/4 else 0);
        end if;
    end for;
    return true, models, covers;
end function;

// Point counts of the fibre product over F_p and F_{p^2} against the trace formula at the first
// three good primes where the model has good reduction.
function trace_formula_agrees(fs, X)
    DN := X`D * X`N;
    checked := 0; p := 3; tried := 0;
    while checked lt 3 and tried lt 60 do
        p := NextPrime(p); tried +:= 1;
        if DN mod p eq 0 then continue; end if;
        Kp := RationalFunctionField(GF(p)); L := Kp;
        try
            for f in fs do
                R<Y> := PolynomialRing(L);
                L := FunctionField(Y^2 - L!Evaluate(PolynomialRing(GF(p))!f, Kp.1));
            end for;
        catch e
            continue;
        end try;
        if Genus(L) ne X`g then continue; end if;
        cnt := [&+[e * #Places(L, e) : e in Divisors(d)] : d in [1..2]];
        if cnt ne [ComputePointsViaTrace(X, p, d) : d in [1..2]] then return false; end if;
        checked +:= 1;
    end while;
    return checked eq 3;
end function;

// A model for X from the file: "hyperelliptic" when a y^2 + h y = f is stored with the genus the
// formula predicts and point counts matching the trace formula; otherwise the double covers of a
// fibre product (stored, or assembled from the covers above W), checked the same way; or "none".
function model_for(X, ModelDir)
    ok, models, covers := stored_models(X`D, X`N, ModelDir);
    if not ok then return "none", _; end if;
    k := Sort(SetToSequence(X`W));
    if IsDefined(models, k) and #models[k] gt 0 then
        e := models[k][1];
        if Type(e[2]) ne MonStgElt then
            f := e[2] + (#e ge 3 select e[3]^2/4 else 0);     // y^2 + h y = f  <=>  (y + h/2)^2 = f + h^2/4
            if Genus(HyperellipticCurve(f)) ne X`g then return "none", _; end if;
            if not trace_formula_agrees([f], X) then return "none", _; end if;
            return "hyperelliptic", _;
        end if;
        fs := CRVDoubleCovers(e[3]);
    else
        ups := AtkinLehnerDoubleCoversOver(X`W, X`D, X`N);
        okg, gens := FibreProductGenerators(X`W, {U : U in ups | IsDefined(covers, U)}, X`D, X`N);
        if not okg then return "none", _; end if;
        failing, classes := AtkinLehnerTripleCheck(covers, X`D, X`N);
        bad := DefinitelyDisagreeingFactors(gens, failing, classes);
        if bad then return "none", _; end if;
        fs := [covers[U] : U in gens];
        if IsDegenerateFactorSet(fs) then return "none", _; end if;
    end if;
    K := FibreProductFunctionField(fs);
    if Genus(K) ne X`g then return "none", _; end if;
    if not trace_formula_agrees(fs, X) then return "none", _; end if;
    return "fibre product", K;
end function;

intrinsic FilterByCanonicalQuadrics(~curves::SeqEnum : ModelDir := "data/models")
{Decides every undecided curve of genus >= 3 for which the canonical embedding is available: at
 D = 1 from the W-fixed weight-2 cusp forms on Gamma_0(N), at D > 1 from a stored model or the
 fibre product of the stored double covers above W (ModelDir).  A stored model y^2 + h y = f is
 hyperelliptic outright.  Decided curves are never overwritten.  Prints one summary line.}
    nq := 0; nm := 0; nhyp := 0; nnon := 0; nskip := 0;
    for i->X in curves do
        if assigned X`IsSubhyp or X`g lt 3 then continue; end if;
        if X`D eq 1 then
            hyp, nrel, g := HyperellipticByQExpansions(X`N, X`W);
            error if g ne X`g, Sprintf("CanonicalQuadrics: X_0(%o)/%o has %o fixed cusp forms but genus %o", X`N, X`W, g, X`g);
            label := Sprintf("CanonicalQuadrics on the W-fixed cusp forms (%o quadrics)", nrel);
            nq +:= 1;
        else
            kind, K := model_for(X, ModelDir);
            if kind eq "none" then nskip +:= 1; continue; end if;
            m0 := X`N eq 1 select "" else "; the model rests on the m = 0 term at N > 1";
            if kind eq "hyperelliptic" then
                hyp := true;
                label := Sprintf("CanonicalQuadrics: a stored hyperelliptic model (genus and point counts checked%o)", m0);
            else
                hyp, nrel := IsHyperellipticFunctionField(K);
                label := Sprintf("CanonicalQuadrics on the stored model (%o quadrics%o)", nrel, m0);
            end if;
            nm +:= 1;
        end if;
        curves[i]`IsSubhyp := hyp;
        curves[i]`IsHyp := hyp;
        curves[i]`TestInWhichProved := label;
        if hyp then nhyp +:= 1; else nnon +:= 1; end if;
    end for;
    printf "CanonicalQuadrics summary: decided %o (hyperelliptic %o, not %o): %o by q-expansions at D = 1, %o by a stored model; %o undecided curves with D > 1 have no usable model\n",
        nhyp + nnon, nhyp, nnon, nq, nm, nskip;
end intrinsic;
