// BiellipticModelCheck.m
//
// Decide between candidate models of a genus-2 bielliptic Atkin-Lehner quotient
// C = X_0(D,N)/W when point counts cannot: two curves with isogenous Jacobians have the SAME
// a_p at every good prime, so the trace formula pins only the isogeny class.  The fixed points
// of the residual Atkin-Lehner involutions, and their fields of definition, are invariants of the
// CURVE together with its AL action (not of its Jacobian), and they do decide.
//
// For an even sextic model  y^2 = A x^6 + B x^4 + C x^2 + D  the group V_4 = <sigma, iota> acts by
//     sigma      : (x,y) -> (-x,  y)   Fix = (0, +-sqrt(D))                 field Q(sqrt(D))
//     sigma*iota : (x,y) -> (-x, -y)   Fix = 2 pts at infty, y/x^3 = +-sqrt(A)  field Q(sqrt(A))
//     iota       : (x,y) -> ( x, -y)   Fix = the 6 Weierstrass points (roots of the sextic)
// and the square classes of A and D are invariants of the model up to x -> lambda x, y -> nu y.
// On the Shimura side the residual group W_full/W acts on C; for each nontrivial coset the number
// of fixed points on C, split by CM discriminant, comes from Ogg's formula (NumFixedPointsByCMOrder,
// CountFixedPointsOnQuotient), and the possible fields of definition of those CM points on C from
// Shimura reciprocity (FieldsOfDefinitionOfCMPoint[Fast], following Gonzalez-Rotger).  A coset
// whose quotient has genus 0 is the hyperelliptic involution iota; a coset with genus-1 quotient is
// one of the two bielliptic involutions, and is matched to sigma or sigma*iota by comparing the a_p
// of its elliptic quotient with the trace formula for X_0(D,N)/<W, coset> -- no hand assignment.
// A candidate is CONSISTENT when, for some a_p-compatible matching, the Galois orbits of the fixed
// points of every AL involution of the model are absorbed by the predicted (disc, count, fields)
// rows (MatchFixedPointOrbits, a multiset check).  The a_p matching is redundant: the fields alone
// decide the pairing (see the note in CheckBiellipticCandidate); it is kept as a cross-check.
//
// Scope: N squarefree (the field-of-definition theory is only implemented there).
//
// Worked out first for X_0(34,3)/w_102 in verify_al_fixed_fields.m; candidate 1 of that entry
// has the right a_p at all primes yet needs fixed points over Q(sqrt(-222)) and Q(sqrt(-74)),
// ramified at 37 -- impossible for CM points of a curve of level 34*3.

declare verbose BiellipticModelCheck, 2;

ALClassRF := recformat< Class : SeqEnum,         // the coset of W, as a sorted sequence of m's
                        Quotient : ShimuraQuot,   // X_0(D,N)/<W, coset>
                        QuotientGenus : RngIntElt,
                        NumFixed : RngIntElt,     // # fixed points of the coset on C
                        Rows : List >;            // [* <disc, count on C, [* possible fields *]> *]

VerdictRF := recformat< D : RngIntElt, N : RngIntElt, Wgens : SetEnum, W : SetEnum,
                        Status : MonStgElt,       // determined | ambiguous | none consistent |
                                                  // not attempted | error
                        Consistent : SeqEnum,     // indices of the consistent candidates
                        IsogenyClass : SeqEnum,   // per candidate: a_p agree with the trace formula
                        Report : MonStgElt >;

// ---------------------------------------------------------------- small helpers

function sqfree(q)
    // squarefree part of a rational number, with sign
    n := Numerator(q)*Denominator(q); s := Sign(n); n := AbsoluteValue(n); r := 1;
    for t in Factorization(n) do if IsOdd(t[2]) then r *:= t[1]; end if; end for;
    return s*r;
end function;

function field_deg(F)
    return Type(F) eq FldRat select 1 else AbsoluteDegree(F);
end function;

function fields_iso(F, G)
    dF := field_deg(F); dG := field_deg(G);
    if dF ne dG then return false; end if;
    if dF eq 1 then return true; end if;
    return IsIsomorphic(AbsoluteField(F), AbsoluteField(G));
end function;

function fldname(F)
    d := field_deg(F);
    if d eq 1 then return "Q"; end if;
    if d eq 2 then return Sprintf("Q(sqrt(%o))", sqfree(Discriminant(MaximalOrder(AbsoluteField(F))))); end if;
    return Sprintf("deg %o, disc %o", d, Factorization(Discriminant(MaximalOrder(AbsoluteField(F)))));
end function;

function quad_orbits(q)
    // Galois orbits of the two points +-sqrt(q): two rational points or one quadratic orbit
    s := sqfree(q);
    return s eq 1 select [* Rationals(), Rationals() *] else [* QuadraticField(s) *];
end function;

// y^2 = A u^3 + B u^2 + C u + D  ==>  Y^2 = U^3 + B U^2 + (A C) U + A^2 D,  U = A u, Y = A y
function ecfromcubic(A, B, C, D)
    return MinimalModel(EllipticCurve([0, B, 0, A*C, A^2*D]));
end function;

function good_primes(f, DN, nprimes)
    // primes of good reduction of the given model y^2 = f, away from DN
    bad := DN * Numerator(Discriminant(f)) * Numerator(LeadingCoefficient(f))
              * LCM([Denominator(c) : c in Coefficients(f)]);
    ps := [];
    p := 3;
    while #ps lt nprimes do
        p := NextPrime(p);
        if bad mod p ne 0 then Append(~ps, p); end if;
    end while;
    return ps;
end function;

// ---------------------------------------------------------------- AL groups and quotients

intrinsic ALGroupFromGenerators(gens::SetEnum, DN::RngIntElt) -> SetEnum
{The subgroup of the Atkin-Lehner group of X_0(D,N), DN = D*N, generated by the w_m, m in gens
 (each m a Hall divisor of DN), as the set of its indices m including 1.}
    for m in gens do
        require DN mod m eq 0 and GCD(m, DN div m) eq 1 : Sprintf("%o is not a Hall divisor of %o", m, DN);
    end for;
    W := {Integers() | 1} join {Integers() | m : m in gens};
    repeat
        n := #W;
        W join:= {AtkinLehnerMul(a, b, DN) : a, b in W};
    until #W eq n;
    return W;
end intrinsic;

intrinsic ALQuotientFromGenerators(D::RngIntElt, N::RngIntElt, gens::SetEnum) -> ShimuraQuot
{The quotient X_0(D,N)/W, W generated by the w_m for m in gens, with its genus assigned.}
    W := ALGroupFromGenerators(gens, D*N);
    X := CreateShimuraQuot(D, N, W);
    X`g := GenusShimuraCurveQuotient(D, N, W);
    return X;
end intrinsic;

// ---------------------------------------------------------------- expected side

intrinsic ExpectedALFixedPointData(X::ShimuraQuot : Fast := true) -> List
{For each nontrivial coset of X`W in the full Atkin-Lehner group of X_0(D,N), the quotient of X by
 it (with genus), the number of its fixed points on X, and those fixed points split by CM
 discriminant together with the possible fields of definition of the CM points on X (Shimura
 reciprocity).  Returns a list of records (Class, Quotient, QuotientGenus, NumFixed, Rows).
 Requires N squarefree.}
    D := X`D; N := X`N; W := X`W; DN := D*N;
    require IsSquarefree(N) : "ExpectedALFixedPointData: the field-of-definition theory needs N squarefree";
    Wfull := {d : d in Divisors(DN) | GCD(d, DN div d) eq 1};

    classes := [];
    seen := {};
    for w in Sort([t : t in Wfull]) do
        if w in seen then continue; end if;
        cl := {AtkinLehnerMul(w, m, DN) : m in W};
        seen join:= cl;
        if cl ne W then Append(~classes, Sort([t : t in cl])); end if;
    end for;

    out := [* *];
    for cl in classes do
        Wq := W join Set(cl);
        Xq := CreateShimuraQuot(D, N, Wq);
        Xq`g := GenusShimuraCurveQuotient(D, N, Wq);
        nfix := Integers()!CountFixedPointsOnQuotient(cl[1], X);
        counts := AssociativeArray();
        for m in cl do
            e := NumFixedPointsByCMOrder(D, N, m);
            for d in Keys(e) do
                counts[d] := (IsDefined(counts, d) select counts[d] else 0) + e[d];
            end for;
        end for;
        rows := [* *];
        total := 0;
        for d in Sort([k : k in Keys(counts)]) do
            if counts[d] eq 0 then continue; end if;
            c := counts[d] / #W;
            require IsIntegral(c) : "fixed-point count on the quotient is not integral";
            c := Integers()!c;
            total +:= c;
            vprintf BiellipticModelCheck, 2 : "  class %o: disc %o, %o fixed point(s) on C; computing fields...\n", cl, d, c;
            fs := Fast select FieldsOfDefinitionOfCMPointFast(X, d) else FieldsOfDefinitionOfCMPoint(X, d);
            Append(~rows, <d, c, fs>);
        end for;
        require total eq nfix : "per-discriminant fixed-point counts do not sum to CountFixedPointsOnQuotient";
        Append(~out, rec< ALClassRF | Class := cl, Quotient := Xq, QuotientGenus := Xq`g,
                                     NumFixed := nfix, Rows := rows >);
    end for;
    return out;
end intrinsic;

// ---------------------------------------------------------------- the comparison

intrinsic MatchFixedPointOrbits(observed::List, expected::List) -> BoolElt
{observed: the Galois orbits of the fixed points of a model involution, one number field (or Q) per
 orbit.  expected: rows <disc, count on the quotient, [* possible fields of definition *]>.  True
 iff the orbits can be assigned to rows so that each row receives exactly its count of points and
 each orbit's field is isomorphic to one of the possible fields of its row.}
    if &+[Integers() | field_deg(F) : F in observed] ne &+[Integers() | r[2] : r in expected] then
        return false;
    end if;
    remaining := [r[2] : r in expected];
    function go(i, rem)
        if i gt #observed then return true; end if;
        F := observed[i]; k := field_deg(F);
        for j in [1..#expected] do
            if rem[j] ge k and &or[fields_iso(F, G) : G in expected[j][3]] then
                rem2 := rem; rem2[j] -:= k;
                if go(i+1, rem2) then return true; end if;
            end if;
        end for;
        return false;
    end function;
    return go(1, remaining);
end intrinsic;

intrinsic CheckBiellipticCandidate(X::ShimuraQuot, f::RngUPolElt, expected::List : NumPrimes := 12) -> BoolElt, MonStgElt, BoolElt
{Test the candidate model y^2 = f (an even sextic) of the genus-2 quotient X against the predicted
 AL fixed-point data (from ExpectedALFixedPointData).  Returns: consistent?, a human-readable
 report, and whether the a_p of the model agree with the trace formula (isogeny class check).}
    D := X`D; N := X`N; DN := D*N;
    require X`g eq 2 : "the quotient must have genus 2";
    require Degree(f) eq 6 and &and[Coefficient(f, i) eq 0 : i in [1, 3, 5]] : "candidate must be an even sextic";
    A := Coefficient(f, 6); B := Coefficient(f, 4); C := Coefficient(f, 2); Dd := Coefficient(f, 0);
    Cv := HyperellipticCurve(f);
    require Genus(Cv) eq 2 : "candidate has the wrong genus";

    rep := Sprintf("y^2 = %o\n", f);
    ps := good_primes(f, DN, NumPrimes);

    // isogeny class: a_p of the model vs the trace formula
    tgt_ap := [p + 1 - ComputePointsViaTrace(X, p, 1) : p in ps];
    obs_ap := [p + 1 - #Points(ChangeRing(Cv, GF(p))) : p in ps];
    isog := tgt_ap eq obs_ap;
    rep cat:= Sprintf("  a_p match the trace formula at p in %o: %o\n", ps, isog);
    if not isog then
        rep cat:= "  *** CONTRADICTION: wrong isogeny class ***\n";
        return false, rep, false;
    end if;

    g0 := [r : r in expected | r`QuotientGenus eq 0];
    g1 := [r : r in expected | r`QuotientGenus eq 1];
    require #g0 + #g1 eq #expected : "an AL involution of a genus-2 curve has quotient of genus 0 or 1";
    require #g0 le 1 : "a genus-2 curve has a unique hyperelliptic involution";
    require #expected le 3 : "the residual AL group of a genus-2 bielliptic curve has order at most 4";

    // model involutions: <label, orbits of fixed points, elliptic quotient or 0>
    model := [* <"sigma      (x,y)->(-x, y)", quad_orbits(Dd), ecfromcubic(A, B, C, Dd)>,
               <"sigma*iota (x,y)->(-x,-y)", quad_orbits(A),  ecfromcubic(Dd, C, B, A)> *];
    iota_orbits := [* Degree(t[1]) eq 1 select Rationals() else NumberField(t[1])
                       : t in Factorization(f) *];

    // a_p-compatible matchings of the genus-1 cosets to the bielliptic model involutions.
    // NB: the traces are NOT needed for the verdict -- the fixed-point fields alone decide the
    // pairing (try both assignments of {sigma, sigma*iota} to the two genus-1 cosets and accept
    // the candidate if either fits, which is what happens below whenever the two elliptic quotients
    // are isogenous).  The a_p matching is kept only as a redundant cross-check that the elliptic
    // quotient attached to a coset lies in the isogeny class of X_0(D,N)/<W, coset>; in the full
    // run over data/bielliptic_candidates.m it rejected nothing.
    q_ap := [[p + 1 - ComputePointsViaTrace(r`Quotient, p, 1) : p in ps] : r in g1];
    m_ap := [[TraceOfFrobenius(mv[3], p) : p in ps] : mv in model];
    matchings := [];   // sequences a with a[i] = index into model for g1[i], injective
    idx := [[j : j in [1..2] | m_ap[j] eq q_ap[i]] : i in [1..#g1]];
    if #g1 eq 0 then
        matchings := [[Integers() | ]];
    elif #g1 eq 1 then
        matchings := [[j] : j in idx[1]];
    else
        matchings := [[j1, j2] : j1 in idx[1], j2 in idx[2] | j1 ne j2];
    end if;
    if #matchings eq 0 then
        rep cat:= Sprintf("  *** CONTRADICTION: no bielliptic involution of the model has the a_p of the genus-1 quotient(s) %o ***\n",
                          [r`Class : r in g1]);
        return false, rep, isog;
    end if;

    function class_name(cl)
        return "{" cat Join([Sprintf("w_%o", m) : m in cl], ",") cat "}";
    end function;
    function rows_name(rows)
        return Join([Sprintf("%o x disc %o: %o", r[2], r[1], [fldname(F) : F in r[3]]) : r in rows], "; ");
    end function;

    consistent := false;
    for a in matchings do
        rep cat:= (#matchings gt 1) select Sprintf("  -- matching %o --\n", a) else "";
        allok := true;
        for i in [1..#g1] do
            mv := model[a[i]]; r := g1[i];
            ok := MatchFixedPointOrbits(mv[2], r`Rows);
            allok and:= ok;
            rep cat:= Sprintf("  %o -> AL coset %o (elliptic quotient, conductor %o)\n", mv[1], class_name(r`Class), Conductor(mv[3]));
            rep cat:= Sprintf("        fixed points observed : %o\n", [fldname(F) : F in mv[2]]);
            rep cat:= Sprintf("        fixed points expected : %o\n", rows_name(r`Rows));
            rep cat:= Sprintf("        %o\n", ok select "MATCH" else "*** CONTRADICTION ***");
        end for;
        if #g0 eq 1 then
            r := g0[1];
            ok := MatchFixedPointOrbits(iota_orbits, r`Rows);
            allok and:= ok;
            rep cat:= Sprintf("  iota       (x,y)->( x,-y) -> AL coset %o (hyperelliptic, quotient genus 0)\n", class_name(r`Class));
            rep cat:= Sprintf("        Weierstrass points    : %o\n", [fldname(F) : F in iota_orbits]);
            rep cat:= Sprintf("        fixed points expected : %o\n", rows_name(r`Rows));
            rep cat:= Sprintf("        %o\n", ok select "MATCH" else "*** CONTRADICTION ***");
        end if;
        consistent or:= allok;
    end for;
    // the third involution of V_4 is not an AL involution when the residual group has order 2
    if #expected eq 1 then
        rep cat:= "  (residual AL group of order 2: only one involution of the model is Atkin-Lehner)\n";
    end if;
    rep cat:= Sprintf("  VERDICT: %o\n", consistent select "CONSISTENT" else "CONTRADICTED");
    return consistent, rep, isog;
end intrinsic;

// ---------------------------------------------------------------- entries of the candidate file

intrinsic ReadBiellipticCandidates(fname::MonStgElt) -> List
{Read a file in the format of fsaia/GenusAtMost2/genus_2_bielliptics_eqn_not_determined.m
 (data/bielliptic_candidates.m): a list of entries [* D, N, W (generating set), [candidates] *].}
    s := Read(fname);
    // the assignment, at the start of a line (the identifier also appears in the header comment)
    i := Position(s, "\ngenus_2_bielliptics_eqn_not_determined");
    require i gt 0 : "not a bielliptic candidate file";
    s := s[i+1..#s];
    j := Position(s, ":=");
    body := s[j+2..#s];
    k := #body;
    while body[k] ne "]" do k -:= 1; end while;
    body := body[1..k];
    return eval ("_<x> := PolynomialRing(Rationals()); return " cat body cat ";");
end intrinsic;

intrinsic CheckBiellipticEntry(e::List : NumPrimes := 12, Fast := true) -> Rec
{Decide an entry [* D, N, W (generating set), [candidate even sextics] *] of the candidate file.
 Returns a record with Status one of "determined" (exactly one consistent candidate), "ambiguous",
 "none consistent", "not attempted" (N not squarefree) or "error", the indices of the consistent
 candidates, per-candidate isogeny-class flags, and a full report.}
    D := e[1]; N := e[2]; gens := e[3]; cands := e[4];
    v := rec< VerdictRF | D := D, N := N, Wgens := gens, Consistent := [], IsogenyClass := [] >;
    if not IsSquarefree(N) then
        v`Status := "not attempted";
        v`Report := Sprintf("X_0(%o,%o)/<%o>: N = %o is not squarefree; field-of-definition theory not implemented\n", D, N, gens, N);
        return v;
    end if;
    try
        X := ALQuotientFromGenerators(D, N, gens);
        v`W := X`W;
        rep := Sprintf("X_0(%o,%o)/W, W = %o (generated by %o), genus %o, %o candidate(s)\n",
                       D, N, Sort([w : w in X`W]), gens, X`g, #cands);
        require X`g eq 2 : "quotient does not have genus 2";
        expected := ExpectedALFixedPointData(X : Fast := Fast);
        for r in expected do
            rep cat:= Sprintf("  AL coset {%o}: %o fixed point(s), quotient genus %o\n",
                              Join([Sprintf("w_%o", m) : m in r`Class], ","), r`NumFixed, r`QuotientGenus);
            for row in r`Rows do
                rep cat:= Sprintf("      %o x disc %o -> %o\n", row[2], row[1], [fldname(F) : F in row[3]]);
            end for;
        end for;
        for i->f in cands do
            vprintf BiellipticModelCheck, 1 : "  candidate %o/%o\n", i, #cands;
            ok, crep, isog := CheckBiellipticCandidate(X, f, expected : NumPrimes := NumPrimes);
            rep cat:= Sprintf("-- candidate %o --\n%o", i, crep);
            Append(~v`IsogenyClass, isog);
            if ok then Append(~v`Consistent, i); end if;
        end for;
        v`Status := #v`Consistent eq 1 select "determined" else
                    (#v`Consistent eq 0 select "none consistent" else "ambiguous");
        rep cat:= Sprintf("==> %o: consistent candidates %o\n", v`Status, v`Consistent);
        v`Report := rep;
    catch err
        v`Status := "error";
        v`Report := Sprintf("X_0(%o,%o)/<%o>: ERROR %o\n", D, N, gens, err`Object);
    end try;
    return v;
end intrinsic;
