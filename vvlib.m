// ---------------------------------------------------------------------------------------------
// Machinery for the vector-valued input form
//     F_f = sum_{gamma in Gamma_0(M)\SL_2(Z)}  (f|_{1/2} gamma) rho(gamma^{-1}) e_0        (GY eq (6))
// evaluated NUMERICALLY.  The point is the CONSTANT term c_0(0): the principal part of F_f only
// sees the cusps where f blows up (oo and 0, GY Lemma 24), but the n = 0 coefficient collects a
// contribution from EVERY coset, so it is not readable off the two scalar q-expansions.
//
// Key design choice that dissolves the metaplectic obstruction recorded in memory
// (phase2-slash-obstruction): we never need a closed-form eta multiplier.  Both the slash and rho
// are evaluated ALONG THE SAME ST-WORD, so they use the same metaplectic lift by construction.
//   S~ = (S, sqrt(tau)),  T~ = (T, 1);   (f|(A,phi))(tau) = phi(tau)^{-1} f(A tau)   (weight 1/2)
//
// WHICH Weil representation: the DUAL of WeilRepresentationST's.  Two independent reasons:
//  (i) well-definedness on cosets at the central element Z = S~^2 = (-I, i).  f|Z = i^{-1} f always
//      (weight 1/2), so the term (f|delta gamma) rho((delta gamma)^{-1}) e_0 is delta-invariant only
//      if rho(Z^{-1}) e_0 = e(1/4) e_0.  WeilRepresentationST gives rho(S)^2 = e(1/4) P (pinned by
//      tests/WeilRepresentation.m), hence rho(Z^{-1}) e_0 = e(-1/4) e_0 -- the wrong sign; the
//      conjugate rep gives e(1/4) e_0.
// (ii) GY Lemma 24 attaches q^{-n/M} to the eta with nrd(eta) in n/M + Z, i.e. exponents = -Q(eta)
//      mod 1, which is rho(T) e_eta = e(-Q(eta)) e_eta -- again the conjugate.
// ---------------------------------------------------------------------------------------------

// --- P^1(Z/M) coset transversal of Gamma_0(M)\SL_2(Z) -----------------------------------------
// Gamma_0(M) g1 = Gamma_0(M) g2  <=>  (c1:d1) = (c2:d2) in P^1(Z/M), so the bottom row classifies.
intrinsic VVCosetReps(M::RngIntElt) -> SeqEnum
{Matrices in SL_2(Z), one per right coset of Gamma_0(M) in SL_2(Z), indexed by P^1(Z/M).}
    seen := {};
    reps := [];
    for c in [0..M-1] do
        for d in [0..M-1] do
            if GCD([c, d, M]) ne 1 then continue; end if;
            // normalize (c:d) up to units mod M
            orb := {@ @};
            for u in [1..M] do
                if GCD(u, M) ne 1 then continue; end if;
                Include(~orb, <(u*c) mod M, (u*d) mod M>);
            end for;
            key := Min([x : x in orb]);
            if key in seen then continue; end if;
            Include(~seen, key);
            // lift (c,d) mod M to a coprime integer pair, then complete the bottom row
            cc := c; dd := d;
            g := GCD(cc, dd);
            if g eq 0 then cc := M; dd := 1; g := 1; end if;
            while GCD(cc, dd) ne 1 do
                // adjust by M to make the pair primitive (always possible since GCD(c,d,M)=1)
                dd +:= M;
            end while;
            _, a, b := XGCD(dd, -cc);      // a*dd - b*cc = 1
            Append(~reps, Matrix(Integers(), 2, 2, [a, b, cc, dd]));
        end for;
    end for;
    return reps;
end intrinsic;

// --- ST-word decomposition ---------------------------------------------------------------------
// Tokens are uniform 2-tuples: <"T", k> = T^k, <"S", 0> = S.  The product, left to right, is g.
intrinsic VVSTWord(g::AlgMatElt) -> SeqEnum
{Write g in SL_2(Z) as a word in S = [0,-1;1,0] and T = [1,1;0,1]; tokens <"T",k> / <"S",0>.}
    require Determinant(g) eq 1 : "Need an element of SL_2(Z).";
    Sm := Matrix(Integers(), 2, 2, [0,-1,1,0]);
    h := g;
    word := [];
    while h[2][1] ne 0 do
        a := h[1][1]; c := h[2][1];
        k := Round(a/c);
        // g = T^k * S * (S^-1 T^-k g)
        Append(~word, <"T", k>);
        Append(~word, <"S", 0>);
        h := Sm^(-1) * Matrix(Integers(), 2, 2, [1,-k,0,1]) * h;
    end while;
    // h = [e, b; 0, e] with e = +-1
    if h[1][1] eq 1 then
        Append(~word, <"T", h[1][2]>);
    else
        Append(~word, <"S", 0>); Append(~word, <"S", 0>);   // -I = S^2
        Append(~word, <"T", -h[1][2]>);
    end if;
    return word;
end intrinsic;

intrinsic VVWordMatrix(word::SeqEnum) -> AlgMatElt
{The SL_2(Z) matrix of an ST-word (product left to right).}
    Sm := Matrix(Integers(), 2, 2, [0,-1,1,0]);
    h := IdentityMatrix(Integers(), 2);
    for t in word do
        if t[1] eq "S" then h := h*Sm; else h := h*Matrix(Integers(),2,2,[1,t[2],0,1]); end if;
    end for;
    return h;
end intrinsic;

// --- numeric evaluation of an eta quotient ------------------------------------------------------
intrinsic VVEtaEval(eta::EtaQuot, z::FldComElt) -> FldComElt
{f(z) for f the eta quotient, evaluated numerically as sum_r c_r prod_d eta(d z)^r_d.}
    R := Parent(eta);
    ds := R`ds;
    CC := Parent(z);
    require Im(z) gt 0 : "z must lie in the upper half plane.";
    etas := [DedekindEta(d*z) : d in ds];
    tot := CC!0;
    for r in Exponents(eta) do
        term := CC!(eta`coeffs[r]);
        for i->d in ds do
            if r[i] ne 0 then term *:= etas[i]^r[i]; end if;
        end for;
        tot +:= term;
    end for;
    return tot;
end intrinsic;

// --- the weight-1/2 metaplectic slash along an ST-word ------------------------------------------
// (f|(A,phi))(tau) = phi(tau)^{-1} f(A tau), with S~ = (S, sqrt(tau)) (principal branch), T~ = (T,1).
// For a word g_1 g_2 ... g_k, f|w = (...((f|g_1)|g_2)...)|g_k, so evaluating at tau walks the word
// BACKWARDS: z_k = tau, z_{i-1} = g_i z_i, and the automorphy factors accumulate along the way.
intrinsic VVSlashEval(eta::EtaQuot, word::SeqEnum, tau::FldComElt) -> FldComElt
{(f|w)(tau) for the weight-1/2 metaplectic slash along the ST-word w.}
    CC := Parent(tau);
    z := tau;
    factor := CC!1;
    for i := #word to 1 by -1 do
        t := word[i];
        if t[1] eq "S" then
            factor /:= Sqrt(z);        // phi_S(z) = sqrt(z), principal branch
            z := -1/z;
        else
            z := z + t[2];
        end if;
    end for;
    return factor * VVEtaEval(eta, z);
end intrinsic;

// --- the Weil representation as complex matrices -------------------------------------------------
// Mirrors WeilRepresentationST but over ComplexField (where 1/sqrt|G| is unambiguously the positive
// real root, so positive_real_sqrt's cyclotomic care is unnecessary).  Dual := true conjugates, i.e.
//   rho*(T) e_eta = e(-nm(eta)) e_eta,   rho*(S)_{i,j} = e(-1/8)/sqrt|G| * e(+<v_i,v_j>)
// which is the representation F_f transforms under (see the header).
intrinsic VVRho(Ld::QuaternionLatticeData, CC::FldCom : Dual := true) -> AlgMatElt, SeqEnum, SeqEnum, RngIntElt
{Complex rho(S) matrix, rho(T) diagonal as a sequence, the ordered discriminant-group elements, and
 the index of the trivial coset.}
    Q := ChangeRing(Ld`Q, Rationals()); dn := Ld`denom;
    elts := [g : g in Ld`disc_grp]; n := #elts;
    vs := [ChangeRing(g@@Ld`to_disc, Rationals()) : g in elts];
    nm := [ (vs[i]*Q, vs[i])/(2*dn^2) : i in [1..n] ];
    i0 := rep{i : i in [1..n] | IsZero(elts[i])};

    ii := CC.1; twopii := 2*Pi(CC)*ii;
    sgn := Dual select -1 else 1;
    ee := func<a | Exp(twopii*CC!(a - Floor(a)))>;

    Tdiag := [ ee(sgn*nm[i]) : i in [1..n] ];
    // Gram of the lifts in one compiled product: Gi[i][j] = (v_i Q, v_j) is an INTEGER (Q integral,
    // the lifts integral), and <v_i,v_j> = Gi[i][j]/dn^2.  So every S-exponent lies in (1/L)Z with
    // L = 2 dn^2, and the whole matrix of root-of-unity indices is one compiled mod-L reduction --
    // which keeps the |G| = 1800 case (3.24e6 entries) down to a table lookup per entry.
    assert &and[&and[IsIntegral(x) : x in Eltseq(v)] : v in vs];
    Vi := Matrix(Integers(), n, 3, &cat[[Integers()!x : x in Eltseq(v)] : v in vs]);
    Gi := Vi * ChangeRing(Q, Integers()) * Transpose(Vi);
    L := 2*dn^2;
    tab := [ Exp(twopii*CC!(k/L)) : k in [0..L-1] ];
    Gidx := ChangeRing((-sgn*2)*Gi, Integers(L));
    c := ee(sgn/8) / Sqrt(CC!n);
    S := ZeroMatrix(CC, n, n);
    for i in [1..n] do
        for j in [1..n] do
            S[i][j] := c * tab[(Integers()!Gidx[i][j]) + 1];
        end for;
    end for;
    return S, Tdiag, elts, i0;
end intrinsic;

// rho(gamma^{-1}) e_0 for gamma given by an ST-word: rho(w) = rho(g_1)...rho(g_k), so
// rho(w^{-1}) e_0 = rho(g_k)^{-1} ... rho(g_1)^{-1} e_0, applied right to left = g_1 first.
// rho is unitary, so rho(S)^{-1} = conj-transpose(S) and rho(T^k)^{-1} = conj of the diagonal.
intrinsic VVRhoInvE0(S::AlgMatElt, Tdiag::SeqEnum, word::SeqEnum, i0::RngIntElt) -> ModTupFldElt
{rho(gamma inverse) applied to e_0, for the ST-word of gamma.}
    CC := BaseRing(S);
    n := Nrows(S);
    v := ZeroMatrix(CC, n, 1); v[i0][1] := 1;
    for t in word do                       // g_1 first: apply rho(g_1)^{-1} first
        if t[1] eq "S" then
            // rho unitary and S symmetric => S^{-1} = conj(S), so S^{-1} v = conj(S conj(v));
            // this avoids storing a second |G| x |G| matrix (415 MB at |G| = 1800, 120 digits).
            for i in [1..n] do v[i][1] := ComplexConjugate(v[i][1]); end for;
            v := S * v;
            for i in [1..n] do v[i][1] := ComplexConjugate(v[i][1]); end for;
        else
            k := t[2];
            for i in [1..n] do v[i][1] *:= (Tdiag[i]^(-k)); end for;
        end if;
    end for;
    return v;
end intrinsic;
