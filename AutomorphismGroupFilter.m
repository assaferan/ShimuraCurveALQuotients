// Automorphism-group (Brandt-Stichtenoth) non-hyperellipticity filter.
//
// Lemma.  Let Y be hyperelliptic of genus g >= 2 with hyperelliptic involution iota, and let
// H <= Aut(Y) with iota notin H.  Then H embeds in Aut(Y)/<iota> <= PGL_2, so
//   * H is cyclic, dihedral (V4 = D_2), A4, S4 or A5;
//   * an elementary abelian H has #H <= 4;
//   * if H ~ V4 then g is odd;
//   * more generally, if H ~ D_n then g = -1 (mod n), or n is odd and g = 0 (mod n).
// References: Brandt-Stichtenoth, manuscripta math. 55 (1986), Prop 2.1 p. 85, Satz 5.1 pp. 89-90,
// Bemerkung 6.3 p. 91; second reference Bujalance-Gamboa-Gromadzki, manuscripta math. 79 (1993),
// Table 2 row 3.a, p. 272.  (The V4 case: lifts of u -> -u and u -> 1/u satisfy
// beta alpha = iota^(g+1) alpha beta.)
//
// Per curve Y = X_0(D,N)/W we compute the known group G_Y: the residual Atkin-Lehner group W_0/W
// together with whichever of S2, V2, V3 descend to Y (they normalize W Gamma_0(DN)), as cosets of
// Q^x W Gamma_0(DN) in the normalizer.  iota is central in Aut(Y), so if iota lies in G_Y it is a
// central involution of G_Y whose quotient has genus 0.  We list those candidates C (for an AL
// involution by GenusShimuraCurveQuotient, for a non-AL one V W_Q by TraceDNewQuotient) and test the
// lemma on every subgroup H <= G_Y that avoids C.  A failure proves Y non-hyperelliptic.
//
// Why the result is trustworthy (prototype: handoff_2026-09-25/scripts/gy.m and the scratchpad
// lemmareach/reach2.m): on the 827 curves undecided after UpdateCurves8 it rules out exactly 9
// (CurveIDs 2555, 2568, 4190, 5635, 5639, 6616, 7926, 7932, 8495, all g = 4), it never rules out
// a recorded-hyperelliptic curve, and the iota-candidate computation finds iota = V2 W8 on
// (1,40)/<w5> and (1,48)/<w3> (quotient genus 0, 2g+2 fixed points).  tests/AutomorphismGroupFilter.m
// pins these values.

// Is g in Q^x * Gamma_0(L)?  (PGL_2: an overall scalar is irrelevant.)
function agInGamma0(g, L)
    ok, c := IsSquare(Determinant(g));
    if not ok then return false; end if;
    h := g/c;
    if not &and[IsIntegral(x) : x in Eltseq(h)] then return false; end if;
    return (Integers()!h[2,1]) mod L eq 0;
end function;

// Is g in Q^x * W * Gamma_0(L)?  Wm is the list of AL matrices of W.
function agInWG(g, Wm, L)
    return exists{m : m in Wm | agInGamma0(g*m^-1, L)};
end function;

// Index of the coset of g among the representatives elts, or 0.
function agFind(g, elts, Wm, L)
    for i->e in elts do
        if agInWG(g*e^-1, Wm, L) then return i; end if;
    end for;
    return 0;
end function;

intrinsic AutomorphismGroupMaxOrder() -> RngIntElt
{Largest order of the known group G_Y that FilterByAutomorphismGroup enumerates; above it the
curve is left undetermined.  The prototype never met a G_Y this large on any curve.}
    return 3000;
end intrinsic;

// The known group G_Y, as a regular permutation group on its elements (cosets of
// Q^x W Gamma_0(DN)).  Also returns a name for each point k of the domain (the element h with
// 1^h = k), and its decomposition: the integer Q if it is the AL coset w_Q, a tuple
// <vname, Q, integer matrix of V> if it is V*w_Q, or 0 if neither was found.
// G is the trivial group if #G_Y exceeds AutomorphismGroupMaxOrder().
function agGroup(X)
    D := X`D; N := X`N; W := X`W; L := D*N;
    M2 := MatrixAlgebra(Rationals(), 2);
    Wm := [M2!al_matrix(w, L) : w in W];
    allQ := [Q : Q in Divisors(L) | GCD(Q, L div Q) eq 1];
    // candidate non-AL generators: <name, matrix over Q, integer matrix for TraceDNewQuotient>
    extra := [];
    if (N mod 4 eq 0) and &and[IsOdd(w) : w in W] then
        S2 := Matrix(Integers(), 2, 2, [2,1,0,2]);
        Append(~extra, <"S2", M2!S2, S2>);
    end if;
    if N mod 8 eq 0 then
        Append(~extra, <"V2", M2!get_V2(L), get_V2(L)>);
    end if;
    if (Valuation(N, 3) eq 2) and ((9 in W) or &and[(w div 3^Valuation(w, 3)) mod 3 eq 1 : w in W]) then
        Append(~extra, <"V3", M2!get_V3(L), get_V3(L)>);
    end if;
    // keep only those that normalize W Gamma_0(L), i.e. descend to Y (a safety net: the conditions
    // above are meant to imply it)
    Wgens := ALsToGens(W, L);
    desc := [];
    for x in extra do
        if &and[agInWG(x[2]*M2!al_matrix(w, L)*x[2]^-1, Wm, L) : w in Wgens] then
            Append(~desc, x);
        else
            vprintf ShimuraQuotients, 1: "KnownAutomorphismGroup: %o does not normalize W on %o\n", x[1], X;
        end if;
    end for;
    gens := [M2!al_matrix(Q, L) : Q in ALsToGens(Seqset(allQ), L)] cat [x[2] : x in desc];
    elts := [M2!1];
    i := 1;
    while i le #elts do
        for x in gens do
            h := elts[i]*x;
            if agFind(h, elts, Wm, L) eq 0 then Append(~elts, h); end if;
        end for;
        i +:= 1;
        if #elts gt AutomorphismGroupMaxOrder() then
            vprintf ShimuraQuotients, 1: "KnownAutomorphismGroup: G_Y too large on %o\n", X;
            return sub<Sym(1) | >, ["1"], [* 1 *];
        end if;
    end while;
    n := #elts;
    Sn := Sym(n);
    perms := [Sn![agFind(elts[j]*x, elts, Wm, L) : j in [1..n]] : x in gens];
    G := sub<Sn | perms>;
    // W Gamma_0(L) is normal in the group generated, so G acts regularly on its cosets
    assert #G eq n;
    named := [<x[1], x[2], x[3]> : x in desc] cat
             [<x[1] cat " " cat y[1], x[2]*y[2], x[3]*y[3]> : x, y in desc | x[1] ne y[1]];
    names := [];
    decomp := [* *];
    for e in elts do
        nm := ""; dc := 0;
        for Q in allQ do
            if agInWG(e*M2!al_matrix(Q, L)^-1, Wm, L) then
                nm := Q eq 1 select "1" else "w" cat IntegerToString(Q); dc := Q; break;
            end if;
        end for;
        if nm eq "" then
            for v in named do
                for Q in allQ do
                    if agInWG(e*(v[2]*M2!al_matrix(Q, L))^-1, Wm, L) then
                        nm := Join(Split(v[1], " "), "*") cat (Q eq 1 select "" else "*w" cat IntegerToString(Q));
                        dc := <v[1], Q, v[3]>;
                        break v;
                    end if;
                end for;
            end for;
        end if;
        if nm eq "" then nm := "?"; end if;
        Append(~names, nm);
        Append(~decomp, dc);
    end for;
    return G, names, decomp;
end function;

// Could the involution with decomposition dc (see agGroup) be the hyperelliptic involution of X?
// True iff its quotient has genus 0, or its genus could not be determined (conservative).
function agMaybeIota(X, dc)
    D := X`D; N := X`N; W := X`W; L := D*N;
    if Type(dc) eq RngIntElt and dc ne 0 then
        WQ := W join {AtkinLehnerMul(dc, w, L) : w in W};
        return GenusShimuraCurveQuotient(D, N, WQ) eq 0;
    elif Type(dc) eq Tup then
        gq := TraceDNewQuotient(dc[3], dc[1], dc[2], W, D, N);
        fix := 2*X`g - 4*gq + 2;
        vprintf ShimuraQuotients, 2: "%o W%o on %o: quotient genus %o, %o fixed points\n", dc[1], dc[2], X, gq, fix;
        // an impossible fixed-point count means we cannot trust gq
        return (gq eq 0) or (fix lt 0) or IsOdd(fix);
    end if;
    return true;    // no decomposition found
end function;

intrinsic KnownAutomorphismGroup(X::ShimuraQuot) -> GrpPerm, SeqEnum
{The known automorphism group G_Y of Y = X: the residual Atkin-Lehner group together with those of
S2, V2, V3 that descend, as a regular permutation group on its elements.  Also returns, for each
point k of the permutation domain, a name for the element h with 1^h = k ("1", "w5", "V3", "V2*w8", ...).
Returns the trivial group if #G_Y exceeds AutomorphismGroupMaxOrder().}
    G, names, _ := agGroup(X);
    return G, names;
end intrinsic;

intrinsic PossibleHyperellipticInvolutions(X::ShimuraQuot) -> SeqEnum
{Names of the central involutions of the known automorphism group of X that could be the
hyperelliptic involution: those whose quotient has genus 0 (or undetermined genus).}
    G, names, decomp := agGroup(X);
    return [names[1^z] : z in Centre(G) | Order(z) eq 2 and agMaybeIota(X, decomp[1^z])];
end intrinsic;

// The subgroups of G that violate the lemma if they avoid iota, smallest first, each with the
// reason.  Conjugacy classes suffice: the conditions are invariant and iota is central.
function agFailingSubgroups(G, g)
    if #G le 2 then return []; end if;
    if IsElementaryAbelian(G) then
        // only V4 (needs g odd) and (Z/2)^3 (never in PGL_2) can fail; larger ones contain a (Z/2)^3
        classes := [H`subgroup : H in Subgroups(G : OrderEqual := 4)];
        if #G ge 8 then classes cat:= [H`subgroup : H in Subgroups(G : OrderEqual := 8)]; end if;
    else
        classes := [H`subgroup : H in Subgroups(G) | H`order ge 4];
    end if;
    fails := [];
    for H in classes do
        if IsCyclic(H) then continue; end if;
        o := #H;
        why := "";
        if o eq 4 then
            if IsEven(g) then why := Sprintf("is V4 and g = %o is even", g); end if;
        elif IsEven(o) and IsIsomorphic(H, DihedralGroup(o div 2)) then
            n := o div 2;
            if not (((g+1) mod n eq 0) or (IsOdd(n) and (g mod n eq 0))) then
                why := Sprintf("is D_%o and g = %o", n, g);
            end if;
        elif (o in {12, 24, 60}) and (IdentifyGroup(H) in {<12,3>, <24,12>, <60,5>}) then
            ;   // A4, S4, A5 embed in PGL_2; no constraint used
        else
            why := Sprintf("is %o, not in PGL_2", IdentifyGroup(H));
        end if;
        if why ne "" then Append(~fails, <H, why>); end if;
    end for;
    Sort(~fails, func<a, b | #a[1] - #b[1]>);
    return fails;
end function;

intrinsic CheckAutomorphismGroup(X::ShimuraQuot) -> BoolElt, MonStgElt
{Returns false if the Brandt-Stichtenoth lemma, applied to the subgroups of the known automorphism
group G_Y that avoid every possible hyperelliptic involution, proves X non-hyperelliptic; in that
case also returns a witness string.  Returns true if nothing is proven.}
    assert X`g ge 2;
    G, names, decomp := agGroup(X);
    fails := agFailingSubgroups(G, X`g);
    if #fails eq 0 then return true, _; end if;
    // Only now decide which central involutions could be iota (the non-AL genus computation is
    // the expensive part), and only for those lying in a failing subgroup.
    Z := Centre(G);
    status := AssociativeArray();
    for f in fails do
        H := f[1];
        avoids := true;
        for z in H do
            if Order(z) ne 2 or z notin Z then continue; end if;
            k := 1^z;
            if not IsDefined(status, k) then status[k] := agMaybeIota(X, decomp[k]); end if;
            if status[k] then avoids := false; break; end if;
        end for;
        if avoids then
            hgens := [names[1^h] : h in Generators(H)];
            return false, Sprintf("H = <%o> %o", Join(hgens, ", "), f[2]);
        end if;
    end for;
    return true, _;
end intrinsic;

intrinsic FilterByAutomorphismGroup(~curves::SeqEnum)
{Mark as non-hyperelliptic the curves whose known automorphism group violates the
Brandt-Stichtenoth lemma (see CheckAutomorphismGroup).}
    for i->X in curves do
        if (i mod 500 eq 0) then
            vprintf ShimuraQuotients, 1: "i = %o/%o\n", i, #curves;
        end if;
        if assigned X`IsSubhyp then continue; end if;
        if X`g lt 3 then continue; end if;
        is_hyp, witness := CheckAutomorphismGroup(X);
        if not is_hyp then
            vprintf ShimuraQuotients, 2: "curve %o: %o\n", X`CurveID, witness;
            curves[i]`IsSubhyp := false;
            curves[i]`IsHyp := false;
            curves[i]`TestInWhichProved := "AutomorphismGroup " cat witness;
        end if;
    end for;
end intrinsic;
