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
// genus formula predicts, AND its point counts over F_p and F_{p^2} at a few good primes have to
// agree with the trace formula.  The second is not redundant: the factors must be written in the
// SAME coordinate on the star, and committed models of one base can come from runs with different
// Hauptmodul normalisations, in which case the compositum is a different curve -- usually of the
// wrong genus, but at X_0(21,2) of the right one, where only the point counts tell.  Before the
// genus, a set of factors also goes through the triple check below, which uses the third double
// cover above the intersection of two of them to decide, where it can, whether two factors are
// written in the same coordinate; a set with two factors that the triples say are not is refused.
// The triples see a coordinate only up to a quadratic twist shared by two covers, so a pass is
// agreement in that sense, and a refusal can be wrong when two covers carry the same wrong twist.
// Inside one pipeline run the first-level stages share a coordinate; the rebase and back-fill
// stages do not (the former substitutes t = (ru+1)/u, the latter stores CurveQuotient's own
// coordinate), so
// AllEquationsAboveCovers passes their keys as Skip and this stage never uses them as factors.

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
    checked := 0; p := 3; tried := 0;
    // run through the good primes until nprimes of them have been checked: a prime where a factor
    // reduces to a square (Y^2 - f reducible) or the genus drops is bad reduction of this model,
    // not a verdict, and is replaced by the next prime rather than counted
    while checked lt nprimes and tried lt 60 do
        p := NextPrime(p); tried +:= 1;
        if DN mod p eq 0 then continue; end if;
        Kp := RationalFunctionField(GF(p));
        L := Kp;
        try
            for f in fs do
                fp := PolynomialRing(GF(p))!f;
                R<Y> := PolynomialRing(L);
                L := FunctionField(Y^2 - L!Evaluate(fp, Kp.1));
            end for;
        catch e
            continue;
        end try;
        if Genus(L) ne X`g then continue; end if;
        cnt := [&+[e * #Places(L, e) : e in Divisors(d)] : d in [1..2]];
        exp := [ComputePointsViaTrace(X, p, d) : d in [1..2]];
        if cnt ne exp then return false, Sprintf("the trace formula disagrees at p = %o", p); end if;
        checked +:= 1;
    end while;
    if checked lt nprimes then return false, Sprintf("only %o good primes among the first 60", checked); end if;
    return true, checked;
end function;

// A set of factors is degenerate when the product of some nonempty subset is a constant times a
// square: if the constant is a square the compositum has the wrong degree (Magma stops with
// "Polynomial must be irreducible"), otherwise the constant field grows and the curve is not
// geometrically irreducible although its genus can still come out right.
intrinsic IsDegenerateFactorSet(fs::SeqEnum[RngUPolElt]) -> BoolElt
{True when the product of some nonempty subset of fs is a constant times a square, i.e. every
 irreducible factor of that product has even multiplicity.  (SquarefreePart returns the radical, so
 it cannot be used for this: SquarefreePart(-3 f^2) is f.)}
    for sub in Subsets({1..#fs}) do
        if IsEmpty(sub) then continue; end if;
        g := &*[fs[i] : i in sub];
        if g eq 0 then return true; end if;
        if forall{fa : fa in Factorization(g) | IsEven(fa[2])} then return true; end if;
    end for;
    return false;
end intrinsic;

// The third index-2 subgroup above U meet V.  Its character is the product of the characters of U
// and V, so it holds the elements that lie in both of them or in neither.
function third_index_two(U, V, full)
    return {Integers()| a : a in full | (a in U) eq (a in V)};
end function;

function class_index(classes, k)
    for i in [1..#classes] do
        if k in classes[i] then return i; end if;
    end for;
    return 0;
end function;

intrinsic AtkinLehnerTripleCheck(covers::Assoc, D::RngIntElt, N::RngIntElt) -> SeqEnum, SeqEnum
{Test whether the stored double covers of one base are written in a single coordinate on the star
 curve.  covers sends an index-2 subgroup U of the full Atkin-Lehner group to the f of its double
 cover y^2 = f(t), with any h^2/4 already added.  If U and V are two of them and T is the third
 index-2 subgroup above U meet V, then f_U f_V f_T is the square of a polynomial whenever the three
 are written in the same coordinate: the three quadratic characters multiply to the trivial one, so
 the product is a square in Q(t), and a square in Q(t) that is a polynomial is the square of one.
 Returns the triples <U, V, T> that fail this, and the partition of the keys into classes joined by
 the triples that pass: two covers in one class are written in the same coordinate up to a
 quadratic twist shared by both, since twisting two of the three by one non-square changes the
 product by its square.}
    DN := D * N;
    full := {Integers()| d : d in Divisors(DN) | GCD(d, DN div d) eq 1};
    ks := [k : k in Keys(covers) | k subset full and #k * 2 eq #full];
    classes := [{k} : k in ks];
    failing := [];
    done := {};
    for i in [1..#ks] do
        for j in [i+1..#ks] do
            T := third_index_two(ks[i], ks[j], full);
            if not IsDefined(covers, T) then continue; end if;
            tri := {ks[i], ks[j], T};
            if tri in done then continue; end if;
            Include(~done, tri);
            g := covers[ks[i]] * covers[ks[j]] * covers[T];
            // every irreducible factor to an even power AND a square leading constant.  The
            // constant matters: rescaling one cover by a non-square rational leaves its zeros where
            // they were and gives a different double cover.  SquarefreePart cannot be used for
            // either part, it returns the radical and so reports -3*f^2 as f.
            sq := g ne 0 and forall{fa : fa in Factorization(g) | IsEven(fa[2])};
            if sq then sq := IsSquare(g); end if;
            if sq then
                a := class_index(classes, ks[i]);
                for k in [ks[j], T] do
                    b := class_index(classes, k);
                    if b ne a then
                        classes[a] join:= classes[b];
                        Remove(~classes, b);
                        a := class_index(classes, ks[i]);
                    end if;
                end for;
            else
                Append(~failing, <ks[i], ks[j], T>);
            end if;
        end for;
    end for;
    return failing, classes;
end intrinsic;

intrinsic DefinitelyDisagreeingFactors(gens::SeqEnum, failing::SeqEnum, classes::SeqEnum)
    -> BoolElt, Tup
{True when two of the double covers gens are known to be written in different coordinates: some
 failing triple holds both of them, and its third member is in the class of one of the two.  That
 one shares a coordinate with the third, so the other cannot, and the two must not be multiplied
 together.  The classes are joined by passing triples, and a triple passes when two of its covers
 carry the same wrong twist; so a correct pair can be refused on that basis.  Against factors
 from runs with different normalisations of t, which is what this is for, the refusal is right.  A failing triple whose third member is in neither class says only that one of the three
 is out of step, not which, so it does not refuse anything.  Also returns the triple.}
    G := Seqset(gens);
    for tri in failing do
        mem := [tri[1], tri[2], tri[3]];
        for i in [1..3] do
            A := mem[(i mod 3) + 1];
            B := mem[((i + 1) mod 3) + 1];
            if A notin G or B notin G then continue; end if;
            cC := class_index(classes, mem[i]);
            if cC ne 0 and (cC eq class_index(classes, A) or cC eq class_index(classes, B)) then
                return true, tri;
            end if;
        end for;
    end for;
    return false, _;
end intrinsic;

intrinsic TripleCheckVerdict(gens::SeqEnum, failing::SeqEnum, classes::SeqEnum) -> MonStgElt
{What the triple check settles about a generating set it did not refuse: whether the generators are
 consistent with one coordinate on the star (a triple twisted by one non-square on two of its covers
 also passes, so this is up to such a shared twist), and whether a failing triple at this base
 involves two of them.}
    G := Seqset(gens);
    extra := "";
    for tri in failing do
        if #({tri[1], tri[2], tri[3]} meet G) ge 2 then
            extra := "; a failing triple at this base involves two of the generators and a cover that is not one of them";
            break;
        end if;
    end for;
    if #gens le 1 then return "conclusive (a single factor, nothing to compare)" cat extra; end if;
    cls := [class_index(classes, g) : g in gens];
    if exists{i : i in [1..#gens] | cls[i] eq 0 or #classes[cls[i]] eq 1} then
        return "inconclusive (a generator is in no triple that passes)" cat extra;
    end if;
    if #Set(cls) eq 1 then
        return "consistent (the generators agree with one coordinate, up to a quadratic twist shared by two of them)" cat extra;
    end if;
    return "inconclusive (the generators fall in classes no triple joins)" cat extra;
end intrinsic;

intrinsic EquationsByFibreProduct(all_eqns::Assoc, all_ws::Assoc, curves::SeqEnum : NPrimes := 3, Skip := {}) -> Assoc, Assoc
{Fill covers that still have no equation by taking the fibre product, over the star curve, of
 their index-2 Atkin-Lehner double covers that do.  A result is kept only when the compositum has
 the genus the Shimura-curve genus formula predicts AND its point counts over NPrimes good primes
 agree with the trace formula.  Reaches covers with no genus-0 quotient, which the Borcherds and
 Schofer route cannot.  The curve stored is the fibre product in its natural affine model, which is
 usually SINGULAR (the double covers share branch points over the base): a model for the function
 field and for point counts through FibreProductFunctionField, not a curve to count or search for
 points on directly.}
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
        cand := [byW[U] : U in ups | IsDefined(byW, U) and byW[U] notin Skip and has_hyperelliptic_eqn(all_eqns, byW[U])];
        if IsEmpty(cand) then continue; end if;
        // factors must share a base: group the candidates by the bases they are written over
        bases := &meet[Keys(all_eqns[j]) : j in cand];
        built := false;
        for b in bases do
            // The triple check sees every stored double cover of the base, not only the ones above
            // W: it is the triples through the other covers that join two factors into one class.
            covers_b := AssociativeArray();
            for j in at_base do
                U := curves[j]`W;
                if #U * 2 ne #full or j in Skip then continue; end if;
                if not IsDefined(all_eqns, j) or not IsDefined(all_eqns[j], b) then continue; end if;
                if Type(all_eqns[j][b]) ne CrvHyp then continue; end if;
                fU, hU := HyperellipticPolynomials(all_eqns[j][b]);
                covers_b[U] := fU + hU^2/4;
            end for;
            failing, classes := AtkinLehnerTripleCheck(covers_b, D, N);
            avail := {curves[j]`W : j in cand | IsDefined(all_eqns[j], b) and Type(all_eqns[j][b]) eq CrvHyp};
            // every generating set is tried, so that a degenerate one does not block the cover
            genlist := [* *];
            k := Ilog2(#full div #X`W);
            if #avail ge k then
                for c in Subsets(avail, k) do
                    I := full; for W2 in c do I := I meet W2; end for;
                    if I eq X`W then Append(~genlist, SetToSequence(c)); end if;
                end for;
            end if;
            for gens in genlist do
            // y^2 + h y = f is (y + h/2)^2 = f + h^2/4: the double cover is y'^2 = f + h^2/4
            fs := [];
            for U in gens do
                f, h := HyperellipticPolynomials(all_eqns[byW[U]][b]);
                Append(~fs, f + h^2/4);
            end for;
            if exists{f : f in fs | f eq 0} then continue; end if;
            if IsDegenerateFactorSet(fs) then
                // correct double covers in one coordinate never form such a set, so this is a
                // defect in the stored equations and is reported like the genus and point-count failures
                printf "WARNING FibreProductCovers: X_0^%o(%o) W=%o over base %o: a subset of the factors multiplies to a constant times a square -- skipped\n",
                    X`D, X`N, Sort(SetToSequence(X`W)), b;
                continue;
            end if;
            bad, tri := DefinitelyDisagreeingFactors(gens, failing, classes);
            if bad then
                printf "WARNING FibreProductCovers: X_0^%o(%o) W=%o over base %o from %o: two of the factors are written in different coordinates (the triple %o fails) -- skipped\n",
                    X`D, X`N, Sort(SetToSequence(X`W)), b, [Sort(SetToSequence(U)) : U in gens],
                    [Sort(SetToSequence(U)) : U in [tri[1], tri[2], tri[3]]];
                continue;
            end if;
            K := FibreProductFunctionField(fs);
            // A rejection is always reported: if the stored double covers are right and share one
            // coordinate, every generating set gives the same curve, so a set that fails here points
            // at a stored double-cover equation that is wrong or in another coordinate.
            if Genus(K) ne X`g then
                printf "WARNING FibreProductCovers: X_0^%o(%o) W=%o over base %o from %o: genus %o, expected %o -- rejected\n",
                    X`D, X`N, Sort(SetToSequence(X`W)), b, [Sort(SetToSequence(U)) : U in gens], Genus(K), X`g;
                continue;
            end if;
            agree, info := trace_formula_agrees(fs, X, NPrimes);
            if not agree then
                printf "WARNING FibreProductCovers: X_0^%o(%o) W=%o over base %o from %o: right genus but %o -- rejected\n",
                    X`D, X`N, Sort(SetToSequence(X`W)), b, [Sort(SetToSequence(U)) : U in gens], info;
                continue;
            end if;
            C := FibreProductCurve(fs);
            if not IsDefined(all_eqns, i) then all_eqns[i] := AssociativeArray(); end if;
            all_eqns[i][b] := C;
            vprintf FibreProductCovers, 1 : "  built W=%o (genus %o) as the fibre product of %o double covers over base %o; trace formula agrees at %o primes; triple check %o\n",
                Sort(SetToSequence(X`W)), X`g, #fs, b, info,
                TripleCheckVerdict(gens, failing, classes);
            built := true;
            break;
            end for;
            if built then break; end if;
        end for;
    end for;
    return all_eqns, all_ws;
end intrinsic;
