// Generalized "complicated AL fixed points on quotient" filter.
//
// Applies prop:complicatedAL (the generalization of Proposition 6 of [FH], which
// TestComplicatedALFixedPointsOnQuotient applies with G an AL group) with G a MIXED group
// containing the non-AL modular involution V2.  This follows Hasegawa's treatment of X_0^*(N):
// the full group of modular involutions replaces W(N) in Prop 6.
//
// Proof being implemented (per curve C = X_0(D,N)/W, genus g_C >= 3, W an AL subgroup with
// W = <W_odd, w_pPart>, p = 2, 8 | N):
//   * G = <W_odd, V2>, V2 = S2 w_pPart S2^-1, is S2 W S2^-1, so C = X/W ~= X/G, and w_N on X/G
//     for N in the pPart-coset of W_odd corresponds to V2 on C.
//   * prop:complicatedAL is applied to G with g1 = w_N1, g2 = w_N2, N1, N2 in that coset.  Its
//     hypotheses (G an elementary 2-group, w_N1, w_N2 notin G, w_N1 w_N2 in G, nu(w_N1) = #G,
//     nu(w_N2) = 3 #G, w_N2 commuting with G and not hyperelliptic on X/G, the congruence and
//     class-number conditions on N2) are checked by CheckComplicatedALHypotheses, which
//     TestComplicatedALFixedPointsOnQuotient uses for the AL case.
//   * V2 must have exactly 4 fixed points on C (the 3 + 1 of the proof); see the docstring.
//
// This is ADDITIVE to FilterByComplicatedALFixedPointsOnQuotient: G contains V2, so it catches
// curves the Atkin-Lehner-only Prop 6 cannot.
//
// It depends on intrinsics from ShimuraQuotients.m, TraceFormula.m and ModularNonALInvolutions.m.

// Level cap: skip curves whose level D*N exceeds this, to avoid runaway trace-formula work
// on a few large highly-composite levels (mirrors NonALModSymMaxLevel for the non-AL filter).
intrinsic GeneralizedComplicatedMaxLevel() -> RngIntElt
{Maximum level D*N for which the generalized complicated-fixed-point filter runs.}
    return 3000;
end intrinsic;

// Session cache for NumFixedPointsNonALOnX: the value depends only on (D, N, vname, Q),
// not on the quotient group W, so it is reused across all curves sharing a level.
GCFP_STORE := NewStore();

function gcfpCache()
    b, A := StoreIsDefined(GCFP_STORE, "cache");
    if not b then A := AssociativeArray(); end if;
    return A;
end function;

// Number of fixed points of a non-AL modular involution u1 = V*W_Q on the FULL curve
// X = X_0(D,N) (before quotienting by any W).  Uses the trace formula:
//   genus(X/(V*W_Q)) = TraceDNewQuotient(V, vname, Q, {1}, D, N),
//   #Fix = 2 g_X - 4 genus(X/(V*W_Q)) + 2   (Riemann-Hurwitz, degree 2).
intrinsic NumFixedPointsNonALOnX(V::AlgMatElt, vname::MonStgElt, Q::RngIntElt,
                                 D::RngIntElt, N::RngIntElt) -> RngIntElt
{Number of fixed points of the non-AL modular involution V*W_Q on X_0(D,N).}
    cache := gcfpCache();
    key := <D, N, vname, Q>;
    cached, val := IsDefined(cache, key);
    if cached then return val; end if;
    gX := GenusShimuraCurve(D, N);
    gQuot := TraceDNewQuotient(V, vname, Q, {Integers() | 1}, D, N);
    val := 2*gX - 4*gQuot + 2;
    cache[key] := val;
    StoreSet(GCFP_STORE, "cache", cache);
    return val;
end intrinsic;

// The single non-AL modular involutions S2, V2, V3 that descend to C = X_0(D,N)/W, returned
// as parallel lists of matrices and names, together with, for each, the set of Atkin-Lehner
// w_Q it does NOT commute with.  Singletons only: the one caller,
// CheckGeneralizedComplicatedFixedPoints, uses nothing else.  (Products of these, and which
// v*W_o are involutions of C, are handled by ModularNonALInvolutionCandidates.)
intrinsic AvailableNonALInvolutions(D::RngIntElt, N::RngIntElt, W::SetEnum)
    -> SeqEnum, SeqEnum, SeqEnum
{Matrices, names, and per-involution non-commuting AL sets for the non-AL modular
involutions on X_0(D,N)/W.}
    DN := D*N;
    Vs := []; V_names := [];
    if (N mod 4 eq 0) and &and[IsOdd(w) : w in W] then
        Append(~Vs, Matrix(Integers(),2,2,[2,1,0,2])); Append(~V_names, "S2");
    end if;
    if (N mod 8 eq 0) then
        Append(~Vs, get_V2(DN)); Append(~V_names, "V2");
    end if;
    if (Valuation(N,3) eq 2) then
        not_commute := false;
        if (9 notin W) then
            not_commute := exists(w){w : w in W | (w div 3^Valuation(w,3)) mod 3 eq 2};
        end if;
        if not not_commute then
            Append(~Vs, get_V3(DN)); Append(~V_names, "V3");
        end if;
    end if;
    all_vs := Vs;
    all_names := V_names;

    als := [Q : Q in Divisors(DN) | GCD(Q, DN div Q) eq 1];
    bad_sets := [];
    for idx->vname in all_names do
        bad := {Integers() |};
        if "S2" in vname then
            bad join:= {w : w in als | IsEven(w)};
        end if;
        // V3 fails to commute with w_m exactly when the 3-free part of m is 2 mod 3 ([FH] Lemma 1;
        // multiplicative in m, as in the descent test above -- e.g. it commutes with w_10).
        if ("V3" in vname) and (9 notin W) then
            bad join:= {w : w in als | (w div 3^Valuation(w,3)) mod 3 eq 2};
        end if;
        Append(~bad_sets, bad);
    end for;
    return all_vs, all_names, bad_sets;
end intrinsic;

// The group checks of prop:complicatedAL.  G = <W_odd, V_p> must equal S_p W S_p^-1: S_p must
// commute with every w in W_odd (modulo Q^* Gamma_0(DN)), and G must be elementary abelian of
// order #W.  Also returns the name of the first failing check, and the matrices of G.
// For V2 (S_p = S2, 8 | N) these checks always pass: S2 = 2T with T = [1,1/2;0,1], T normalises
// Gamma_0(L) when 4 | L and takes each odd w_m to a matrix of the same AL form, so G = T W T^-1.
// They reject for V3 (e.g. (10,153)/<2,9,85>, where S3 does not commute with w2 or w85).
intrinsic GeneralizedComplicatedMixedGroup(DN::RngIntElt, Wodd::SetEnum, V::AlgMatElt,
                                           Sp::AlgMatElt) -> BoolElt, MonStgElt, SeqEnum
{True iff G = W_odd cat V*W_odd is S_p W S_p^-1 modulo Q^* Gamma_0(DN) and elementary abelian of
order 2 #W_odd; otherwise false and the failing check ("commute", "square", "abelian" or
"distinct").  The third value is the sequence of matrices of G.}
    M2Q := MatrixAlgebra(Rationals(), 2);
    Sp := M2Q!Sp;
    Vq := M2Q!V;
    Wmats := [M2Q!al_matrix(w, DN) : w in Wodd];
    Gmats := Wmats cat [Vq*A : A in Wmats];
    if not &and[IsTrivialModGamma0(Sp*A*Sp^-1*A^-1, DN) : A in Wmats] then
        return false, "commute", Gmats;
    end if;
    if not &and[IsTrivialModGamma0(x*x, DN) : x in Gmats] then return false, "square", Gmats; end if;
    if not &and[IsTrivialModGamma0(x*y*x^-1*y^-1, DN) : x, y in Gmats] then
        return false, "abelian", Gmats;
    end if;
    if exists{i : i, j in [1..#Gmats] | i lt j and IsTrivialModGamma0(Gmats[i]*Gmats[j]^-1, DN)} then
        return false, "distinct", Gmats;
    end if;
    return true, "", Gmats;
end intrinsic;

intrinsic CheckGeneralizedComplicatedFixedPoints(X::ShimuraQuot) -> BoolElt, MonStgElt
{Returns true and a witness string if C = X_0(D,N)/W is proven non-hyperelliptic by
prop:complicatedAL applied to a MIXED group G generated by W_odd (the w in W coprime to p)
together with the non-AL modular involution V_p = S_p w_pPart S_p^-1 at p = 2 (V2,
S2 = [2,1;0,2]; V3 is excluded, see below).  Via the isomorphism C = X/W = X/G the involution
w_pPart on X/G plays the role of g2 = w_N2 and g1 = w_N1, for AL involutions N1, N2 lying OUTSIDE
G, so this reaches the star/full-W quotients the pure-AL test cannot.  N2 supplies three
conjugate complicated CM fixed points and N1 the rational fourth.  The hypotheses of the
proposition are checked by CheckComplicatedALHypotheses; genus(C) >= 3 is checked here.
SOUNDNESS: in addition, V_p must have exactly 4 fixed points on C (the 3 + 1 of the proof).
On a hyperelliptic curve every involution other than the hyperelliptic involution has 0/2/4
fixed points, and the hyperelliptic involution itself has 2g+2; so if the true number of fixed
points of V_p on C equals 2g+2 (e.g. X_0(35,16)) then V_p is itself hyperelliptic and C is
hyperelliptic, and any value other than 4 means the 3 + 1 count is inflated by non-AL coset
contributions and the argument breaks.  We therefore require the true count, namely (1/#W)
times the sum over w in W of NumFixedPoints of V_p*w on X, to equal 4.}
    if X`g lt 3 then return false, _; end if;
    D := X`D; N := X`N; W := X`W; DN := D*N;
    // Non-AL modular involutions exist only when 4 | N or 9 || N.
    if (N mod 4 ne 0) and (Valuation(N,3) ne 2) then return false, _; end if;
    if DN gt GeneralizedComplicatedMaxLevel() then return false, _; end if;

    als := [Q : Q in Divisors(DN) | GCD(Q, DN div Q) eq 1];
    all_vs, all_names, _ := AvailableNonALInvolutions(D, N, W);

    for idx->vname in all_names do
        // V2 only.  V3 is excluded because a V3 certificate can never be valid: Prop 6 needs
        // G = <W_odd, V3> = S3 W S3^-1 elementary abelian with w_N1, w_N2 notin G, which forces
        // S3 to commute with every w in W_odd; but then Fix(V3 on C) >= 8, contradicting the
        // guard Fix(V3 on C) = 4 below.  (When S3 does not commute, V3 w_m V3^-1 = w_{9m} for
        // m = 2 mod 3, so G contains w_9, is non-abelian of order 16, and N1, N2 lie in G.)
        if vname ne "V2" then continue; end if;
        p := 2;
        V := all_vs[idx];

        pPart := p^Valuation(DN, p);                    // the AL involution V_p replaces
        if pPart notin W then continue; end if;         // V_p replaces w_pPart, so it must be in W
        Wodd := {w : w in W | GCD(w, p) eq 1};          // Atkin-Lehner part of G
        if #W ne 2*#Wodd then continue; end if;         // W = <W_odd, w_pPart>, so X/W ~= X/G
        M2Q := MatrixAlgebra(Rationals(), 2);
        ok := GeneralizedComplicatedMixedGroup(DN, Wodd, V, M2Q![p, 1, 0, p]);
        if not ok then continue; end if;

        // SOUNDNESS GUARD: V_p must have exactly 4 fixed points on C (not 2g+2, not inflated).
        nuVp := (&+[NumFixedPointsNonALOnX(V, vname, w, D, N) : w in W]) / #W;
        if nuVp ne 4 then continue; end if;

        // N2: three conjugate complicated fixed points; N1: the rational fourth.  Both AL
        // involutions in the pPart-coset of W_odd (i.e. outside G but representing w_pPart).
        // The hypotheses of prop:complicatedAL are checked by CheckComplicatedALHypotheses.
        for N2 in als do
            if N2 eq 1 or N2 in Wodd then continue; end if;
            if AtkinLehnerMul(N2, pPart, DN) notin Wodd then continue; end if;
            for N1 in als do
                if N1 eq 1 or N1 in Wodd or N1 eq N2 then continue; end if;
                if AtkinLehnerMul(N1, pPart, DN) notin Wodd then continue; end if;
                if CheckComplicatedALHypotheses(D, N, Wodd, N1, N2 : V := V) then
                    return true, Sprintf(
                        "GeneralizedComplicatedFixedPoints: %o, pPart=%o, N1=%o (nu=%o), N2=%o (nu=%o)",
                        vname, pPart, N1, #W, N2, 3*#W);
                end if;
            end for;
        end for;
    end for;
    return false, _;
end intrinsic;

intrinsic FilterByGeneralizedComplicatedFixedPoints(~curves::SeqEnum)
{Mark curves proven non-hyperelliptic by the generalized Prop 6 (non-AL modular involution
as the free-action certifier). Additive to FilterByComplicatedALFixedPointsOnQuotient.}
    for i->X in curves do
        if assigned X`IsSubhyp then continue; end if;
        if X`g lt 3 then continue; end if;
        ok, witness := CheckGeneralizedComplicatedFixedPoints(X);
        if ok then
            curves[i]`IsSubhyp := false;
            curves[i]`IsHyp := false;
            curves[i]`TestInWhichProved := witness;
        end if;
        if (i mod 200 eq 0) then
            vprintf ShimuraQuotients, 1: "i = %o/%o\n", i, #curves;
        end if;
    end for;
end intrinsic;
