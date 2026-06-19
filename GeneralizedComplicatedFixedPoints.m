// Generalized "complicated AL fixed points on quotient" filter.
//
// Upgrades Proposition 6 of [FH] (implemented in TestComplicatedALFixedPointsOnQuotient)
// by letting the free-action certifier -- the "N1" of the [FH] proof -- be a NON-AL
// modular involution in the full modular involution group G, instead of only an
// Atkin-Lehner involution.  This follows Hasegawa's treatment of X_0^*(N): when
// X_0(D,N)/W carries extra non-AL modular involutions (S2, V2, V3 and products),
// the full group G of modular involutions replaces W(N) in Prop 6.
//
// Proof being implemented (per curve C = X_0(D,N)/W, genus g_C >= 3, W an AL subgroup):
//   * Pick an Atkin-Lehner involution w_{N2} with N2 | D*N, gcd(N2, D*N/N2)=1, N2 notin W,
//     h(-4 N2) = 0 (mod 3) and the [FH] congruence conditions.  These force Fix(w_{N2})
//     on X to consist solely of CM points of the single order of discriminant -4 N2, of
//     size nu(w_{N2}) = NumFixedPoints(D,N,N2), with prime divisors in {2,3} and divisible
//     by 3.  We further require nu(w_{N2}) = 3*#W exactly, so that the elementary 2-group W
//     (order #W) acting freely on Fix(w_{N2}) yields exactly 3 fixed points of w_{N2} on C.
//   * The free action is certified by an involution u1 in G that commutes with w_{N2} and
//     has nu(u1) on X = #W = nu(w_{N2})/3.  In the [FH] (Atkin-Lehner) setting u1 is the
//     AL involution N1 with N1*W = N2*W; here u1 is allowed to be a non-AL modular
//     involution V*W_Q.  (The disjointness of the fixed sets -- i.e. freeness -- follows
//     from the CM structure of Fix(w_{N2}).)
//   * Three fixed points is odd, so by Ogg's criterion C is not hyperelliptic.
//
// This is ADDITIVE to FilterByComplicatedALFixedPointsOnQuotient: it only looks for a
// NON-AL partner u1, catching curves the Atkin-Lehner-only Prop 6 cannot.
//
// This file is a STANDALONE package (deliberately NOT in ShimuraQuotients.spec): Attach it
// on top of the spec.  It depends on intrinsics from ShimuraQuotients.m, TraceFormula.m and
// ModularNonALInvolutions.m.  Once validated it can be folded into ModularNonALInvolutions.m
// and added to the spec.

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

// The non-AL modular involutions available on C = X_0(D,N)/W, returned as parallel lists
// of matrices and names, together with, for each, the set of Atkin-Lehner w_Q it does NOT
// commute with (so that V*W_Q is a well-defined involution and V commutes with w_Q).
// Mirrors the Vs/bad_ws logic of CheckModularNonALInvolutionTrace.
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
    all_vs := Vs cat [v1*v2 : v1, v2 in Vs | v1 ne v2];
    all_names := V_names cat [v1 cat " " cat v2 : v1, v2 in V_names | v1 ne v2];

    als := [Q : Q in Divisors(DN) | GCD(Q, DN div Q) eq 1];
    bad_sets := [];
    for idx->vname in all_names do
        bad := {Integers() |};
        if "S2" in vname then
            bad join:= {w : w in als | IsEven(w)};
        end if;
        if ("V3" in vname) and (9 notin W) then
            bad join:= {w : w in als | exists(p){p : p in PrimeDivisors(w) |
                                                 (p^Valuation(w,p) mod 3) eq 2}};
        end if;
        Append(~bad_sets, bad);
    end for;
    return all_vs, all_names, bad_sets;
end intrinsic;

// Atkin-Lehner N2 candidates carrying the [FH] CM conditions, with nu(w_{N2}) = 3*#W.
intrinsic GeneralizedN2Candidates(D::RngIntElt, N::RngIntElt, W::SetEnum) -> SeqEnum
{Atkin-Lehner divisors N2 of D*N satisfying the [FH] Prop 6 CM conditions with
NumFixedPoints(D,N,N2) = 3*#W, and N2 notin W.}
    DN := D*N;
    sizeW := #W;
    N2s := [];
    for N2 in Divisors(DN) do
        if N2 eq 1 then continue; end if;
        if GCD(N2, DN div N2) ne 1 then continue; end if;
        if N2 in W then continue; end if;
        if ClassNumber(-4*N2) mod 3 ne 0 then continue; end if;
        if not ((N2 mod 4 ne 3) or ((N2 mod 8 eq 3) and IsEven(N))
                                 or ((N2 mod 8 eq 7) and IsEven(D))) then
            continue;
        end if;
        nu2 := NumFixedPoints(D, N, N2);
        if nu2 eq 0 then continue; end if;
        if nu2 mod 3 ne 0 then continue; end if;
        if not (PrimeDivisors(nu2) subset [2,3]) then continue; end if;
        if nu2 ne 3*sizeW then continue; end if;   // exactly 3 fixed points on C
        Append(~N2s, N2);
    end for;
    return N2s;
end intrinsic;

intrinsic CheckGeneralizedComplicatedFixedPoints(X::ShimuraQuot) -> BoolElt, MonStgElt
{Returns true and a witness string if C = X_0(D,N)/W is proven non-hyperelliptic by the
generalized [FH] Prop 6.  The Atkin-Lehner group is replaced by a MIXED group G generated by
W_odd (the w in W coprime to p) together with a non-AL modular involution V_p at a split prime
p dividing N.  Via the isomorphism C = X/W = X/G the involution w_pPart on X/G plays the Prop 6
role with a pair of Atkin-Lehner involutions N1, N2 lying OUTSIDE G, so this reaches the
star/full-W quotients the pure-AL test cannot.  N2 supplies three conjugate complicated CM
fixed points and N1 the rational fourth.
SOUNDNESS: the certificate is valid only when V_p has exactly 4 fixed points on C (the 3 + 1
of the proof).  On a hyperelliptic curve every involution other than the hyperelliptic
involution has 0/2/4 fixed points, and the hyperelliptic involution itself has 2g+2; so if the
true number of fixed points of V_p on C equals 2g+2 (e.g. X_0(35,16)) then V_p is itself
hyperelliptic and C is hyperelliptic, and any value other than 4 means the 3 + 1 count is
inflated by non-AL coset contributions and the argument breaks.  We therefore require the true
count, namely (1/#W) times the sum over w in W of NumFixedPoints of V_p*w on X, to equal 4.}
    if X`g lt 3 then return false, _; end if;
    D := X`D; N := X`N; W := X`W; DN := D*N;
    // Non-AL modular involutions exist only when 4 | N or 9 || N.
    if (N mod 4 ne 0) and (Valuation(N,3) ne 2) then return false, _; end if;
    if DN gt GeneralizedComplicatedMaxLevel() then return false, _; end if;

    als := [Q : Q in Divisors(DN) | GCD(Q, DN div Q) eq 1];
    all_vs, all_names, bad_sets := AvailableNonALInvolutions(D, N, W);

    for idx->vname in all_names do
        if #vname gt 2 then continue; end if;          // single involutions V2/V3 only
        if   "2" in vname then p := 2;
        elif "3" in vname then p := 3;
        else continue; end if;
        V := all_vs[idx]; bad := bad_sets[idx];

        pPart := p^Valuation(DN, p);                    // the AL involution V_p replaces
        if pPart notin W then continue; end if;         // V_p replaces w_pPart, so it must be in W
        Wodd := {w : w in W | GCD(w, p) eq 1};          // Atkin-Lehner part of G
        if #W ne 2*#Wodd then continue; end if;         // W = <W_odd, w_pPart>, so X/W ~= X/G
        if #(Wodd meet bad) ne 0 then continue; end if; // V_p must commute with all of W_odd

        // SOUNDNESS GUARD: V_p must have exactly 4 fixed points on C (not 2g+2, not inflated).
        nuVp := (&+[NumFixedPointsNonALOnX(V, vname, w, D, N) : w in W]) / #W;
        if nuVp ne 4 then continue; end if;

        // N2: three conjugate complicated fixed points; N1: the rational fourth.  Both AL
        // involutions in the pPart-coset of W_odd (i.e. outside G but representing w_pPart).
        for N2 in als do
            if N2 eq 1 or N2 in Wodd then continue; end if;
            if AtkinLehnerMul(N2, pPart, DN) notin Wodd then continue; end if;
            if ClassNumber(-4*N2) mod 3 ne 0 then continue; end if;
            if not ((N2 mod 4 ne 3) or ((N2 mod 8 eq 3) and IsEven(N))
                                     or ((N2 mod 8 eq 7) and IsEven(D))) then continue; end if;
            nfixed := NumFixedPoints(D, N, N2);
            if nfixed eq 0 or (nfixed mod 3 ne 0) or not (PrimeDivisors(nfixed) subset [2,3]) then continue; end if;
            target := 2^Valuation(nfixed, 2);
            if 2*#Wodd ne target then continue; end if;
            for N1 in als do
                if N1 eq 1 or N1 in Wodd or N1 eq N2 then continue; end if;
                if AtkinLehnerMul(N1, pPart, DN) notin Wodd then continue; end if;
                if NumFixedPoints(D, N, N1) ne target then continue; end if;
                if AtkinLehnerMul(N1, N2, DN) in Wodd then
                    return true, Sprintf(
                        "GeneralizedComplicatedFixedPoints: %o, pPart=%o, N1=%o (nu=%o), N2=%o (nu=%o)",
                        vname, pPart, N1, target, N2, nfixed);
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
