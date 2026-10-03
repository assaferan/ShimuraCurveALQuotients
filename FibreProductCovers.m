// FibreProductCovers.m
//
// A cover Y = X_0(D,N)/W is Galois over the star curve, with group W_full/W an elementary
// abelian 2-group.  Its function field is therefore the compositum of its quadratic
// subextensions, and those correspond to the index-2 subgroups of W_full containing W, that is
// to the covers X_0(D,N)/W'' with [W_full : W''] = 2.  Each of those is a double cover of the
// star, so it carries an equation y^2 = f(t) in the star's coordinate.  Writing k for
// log_2 [W_full : W], any k of them whose groups meet in exactly W generate the function field of
// Y, so Y is their fibre product over the star:
//
//     y_1^2 = f_1(t),  ...,  y_k^2 = f_k(t).
//
// What this buys: the construction never needs a quotient of genus 0, so it reaches curves that
// are not subhyperelliptic, which is exactly where the Borcherds and Schofer route stops.  It
// builds the top curve X_0(D,N) itself whenever enough of the big Atkin-Lehner quotients are
// known.
//
// What accepts a result: the genus of the compositum has to equal the genus the Shimura-curve
// genus formula predicts.  That is not automatic.  The factors must be written in the SAME
// coordinate on the star, and committed models of one base can come from runs with different
// Hauptmodul normalisations, in which case the compositum is a different curve and its genus
// comes out wrong.  A caller that pools equations from several runs must treat the genus as the
// test it is, not as a formality.

declare verbose FibreProductCovers, 2;

// y^2 = f(x) with deg f = d is a double cover of the x-line; in the weighted ambient the fibre
// coordinate has weight ceil(d/2), and the right-hand side homogenises to degree twice that.
function fibre_weight(f)
    return Ceiling(Degree(f) / 2);
end function;

function homogenise(f, w, s, z)
    return &+[Coefficient(f, i) * s^i * z^(2*w - i) : i in [0..Degree(f)]];
end function;

intrinsic FibreProductFunctionField(fs::SeqEnum) -> FldFun
{The function field of the fibre product over the t-line of the double covers y_i^2 = f_i(t).}
    require not IsEmpty(fs) : "need at least one factor";
    FF := RationalFunctionField(Rationals());
    t := FF.1;
    K := FF;
    for f in fs do
        R<Y> := PolynomialRing(K);
        K := FunctionField(Y^2 - K!Evaluate(f, t));
    end for;
    return K;
end intrinsic;

intrinsic FibreProductCurve(fs::SeqEnum) -> Crv, SeqEnum
{The fibre product over the t-line of the double covers y_i^2 = f_i(t), as a curve in a weighted
 projective space with coordinates s, z of weight 1 and one fibre coordinate per factor.  Also
 returns its defining polynomials.}
    require not IsEmpty(fs) : "need at least one factor";
    ws := [fibre_weight(f) : f in fs];
    names := ["s", "z"] cat ["y" cat IntegerToString(i) : i in [1..#fs]];
    Pamb := WeightedProjectiveSpace(Rationals(), [1, 1] cat ws);
    AssignNames(~Pamb, names);
    cs := [Pamb.i : i in [1..#names]];
    s := cs[1]; z := cs[2];
    eqns := [cs[2+i]^2 - homogenise(fs[i], ws[i], s, z) : i in [1..#fs]];
    return Curve(Pamb, eqns), eqns;
end intrinsic;

// The index-2 subgroups of the full Atkin-Lehner group that contain W.  Each one's quotient is a
// double cover of the star curve.
function index_two_over(W, full, DN)
    out := {};
    for S in Subsets(full diff W) do
        if IsEmpty(S) then continue; end if;
        Wd := {Integers()| AtkinLehnerMul(a, b, DN) : a in W join S, b in W join S};
        if #Wd eq (#full div 2) then Include(~out, Wd); end if;
    end for;
    return out;
end function;

function has_hyperelliptic_eqn(all_eqns, i)
    return IsDefined(all_eqns, i) and exists{b : b in Keys(all_eqns[i]) | Type(all_eqns[i][b]) eq CrvHyp};
end function;

intrinsic AtkinLehnerDoubleCoversOver(W::SetEnum, D::RngIntElt, N::RngIntElt) -> SetEnum
{The Atkin-Lehner subgroups of index 2 in the full group that contain W.  The quotient of
 X_0(D,N) by each of them is a double cover of the star curve, and X_0(D,N)/W is the fibre
 product of enough of them.}
    DN := D * N;
    full := {Integers()| d : d in Divisors(DN) | GCD(d, DN div d) eq 1};
    require W subset full : "W must be a subgroup of the Atkin-Lehner group";
    return index_two_over(W, full, DN);
end intrinsic;

intrinsic FibreProductGenerators(W::SetEnum, avail::SetEnum, D::RngIntElt, N::RngIntElt)
    -> BoolElt, SeqEnum
{Choose index-2 subgroups from avail whose intersection is exactly W, so that the corresponding
 double covers of the star generate the function field of X_0(D,N)/W.  Returns false when avail
 does not contain such a set.}
    DN := D * N;
    full := {Integers()| d : d in Divisors(DN) | GCD(d, DN div d) eq 1};
    k := Ilog2(#full div #W);
    if #avail lt k then return false, _; end if;
    for c in Subsets(avail, k) do
        I := full;
        for W2 in c do I := I meet W2; end for;
        if I eq W then return true, SetToSequence(c); end if;
    end for;
    return false, _;
end intrinsic;

// Point counts of the compositum over F_p and F_{p^2} against the Eichler-Selberg trace formula on
// the W-fixed part of the D-new space.  The genus check cannot tell a fibre product of factors in
// different coordinates from the right curve when the degrees happen to agree; this can.
function trace_formula_agrees(fs, X, nprimes)
    DN := X`D * X`N;
    ps := [p : p in [5, 7, 11, 13, 17, 19, 23, 29] | DN mod p ne 0][1..nprimes];
    checked := 0;
    for p in ps do
        Kp := RationalFunctionField(GF(p));
        L := Kp;
        // A factor can reduce to a square mod p, making Y^2 - f reducible; that is bad reduction
        // of this model at p, not a verdict, so the prime is skipped rather than raised.
        try
            for f in fs do
                fp := PolynomialRing(GF(p))!f;
                R<Y> := PolynomialRing(L);
                L := FunctionField(Y^2 - L!Evaluate(fp, Kp.1));
            end for;
        catch e
            continue;
        end try;
        if Genus(L) ne X`g then continue; end if;          // bad reduction of this model
        cnt := [&+[e * #Places(L, e) : e in Divisors(d)] : d in [1..2]];
        exp := [ComputePointsViaTrace(X, p, d) : d in [1..2]];
        if cnt ne exp then return false, p; end if;
        checked +:= 1;
    end for;
    return checked gt 0, checked;
end function;

intrinsic EquationsByFibreProduct(all_eqns::Assoc, all_ws::Assoc, curves::SeqEnum : NPrimes := 3) -> Assoc, Assoc
{Fill covers that still have no equation by taking the fibre product, over the star curve, of
 their index-2 Atkin-Lehner double covers that do.  A result is kept only when the compositum has
 the genus the Shimura-curve genus formula predicts AND its point counts over NPrimes good primes
 agree with the trace formula.  Reaches covers with no genus-0 quotient, which the Borcherds and
 Schofer route cannot.}
    labels := [k : k in Keys(all_eqns)];
    if IsEmpty(labels) then return all_eqns, all_ws; end if;
    X0 := curves[labels[1]];
    D := X0`D; N := X0`N; DN := D*N;
    full := {Integers()| d : d in Divisors(DN) | GCD(d, DN div d) eq 1};
    // every cover at this base, by W
    at_base := [i : i in [1..#curves] | curves[i]`D eq D and curves[i]`N eq N];
    byW := AssociativeArray();
    for i in at_base do byW[curves[i]`W] := i; end for;
    for i in at_base do
        X := curves[i];
        if X`W eq full then continue; end if;
        if IsDefined(all_eqns, i) and not IsEmpty(Keys(all_eqns[i])) then continue; end if;
        ups := index_two_over(X`W, full, DN);
        cand := [byW[U] : U in ups | IsDefined(byW, U) and has_hyperelliptic_eqn(all_eqns, byW[U])];
        if IsEmpty(cand) then continue; end if;
        // factors must share a base: group the candidates by the bases they are written over
        bases := &meet[Keys(all_eqns[j]) : j in cand];
        built := false;
        for b in bases do
            avail := {curves[j]`W : j in cand | IsDefined(all_eqns[j], b) and Type(all_eqns[j][b]) eq CrvHyp};
            ok, gens := FibreProductGenerators(X`W, avail, D, N);
            if not ok then continue; end if;
            fs := [HyperellipticPolynomials(all_eqns[byW[U]][b]) : U in gens];
            if exists{f : f in fs | f eq 0} then continue; end if;
            K := FibreProductFunctionField(fs);
            if Genus(K) ne X`g then
                vprintf FibreProductCovers, 1 : "  fibre product for W=%o over base %o has genus %o, expected %o: rejected\n",
                    Sort(SetToSequence(X`W)), b, Genus(K), X`g;
                continue;
            end if;
            agree, info := trace_formula_agrees(fs, X, NPrimes);
            if not agree then
                vprintf FibreProductCovers, 1 : "  fibre product for W=%o over base %o has the right genus but fails the trace formula at p=%o: rejected\n",
                    Sort(SetToSequence(X`W)), b, info;
                continue;
            end if;
            C := FibreProductCurve(fs);
            if not IsDefined(all_eqns, i) then all_eqns[i] := AssociativeArray(); end if;
            all_eqns[i][b] := C;
            vprintf FibreProductCovers, 1 : "  built W=%o (genus %o) as the fibre product of %o double covers over base %o; trace formula agrees at %o primes\n",
                Sort(SetToSequence(X`W)), X`g, #fs, b, info;
            built := true;
            break;
        end for;
    end for;
    return all_eqns, all_ws;
end intrinsic;
