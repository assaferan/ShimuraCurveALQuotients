// The vector-valued Borcherds input form
//     F_f = sum_{gamma in Gamma_0(M)\SL_2(Z)}  (f |_{1/2} gamma)  rho(gamma^{-1}) e_0        [GY, eq (6)]
// evaluated numerically, and the constant terms c_eta(0) that Schofer's m = 0 term is built from.
//
// WHY THIS EXISTS.  The PRINCIPAL part of F_f is given in closed form by [GY, Lemma 24] in terms of
// f's two scalar q-expansions at the cusps oo and 0 -- that is what SchoferFormula.m consumes, and it
// is validated by the m > 0 machinery.  The CONSTANT term is not: the principal part only sees the
// cusps where f blows up, whereas the n = 0 coefficient of the coset sum collects a contribution from
// EVERY coset.  Measured here, it turns out that
//     c_0(0) = 0        and       the m = 0 multiplier = (1/2) c_eta(0) at any NONZERO ISOTROPIC eta
// (all such eta carry the same value).  That is why no linear combination of Coefficient(foo,0) and
// Coefficient(f0,0) ever reproduced the measured multipliers.
//
// THE METAPLECTIC POINT.  Getting f|gamma with the correct metaplectic phase looks like it needs an
// explicit generalized-eta multiplier system.  It does not, as long as one never asks for a closed
// form: decompose each coset representative into an ST-WORD and evaluate BOTH the weight-1/2 slash
// and rho along the SAME word.  They then share one metaplectic lift by construction, so the +-1
// ambiguity of the double cover never arises.  Conventions:
//     S~ = (S, sqrt(tau))  (principal branch),   T~ = (T, 1),
//     (f | (A,phi))(tau) = phi(tau)^{-1} f(A tau)                        (weight 1/2)
//
// WHICH Weil representation: the DUAL of WeilRepresentationST's.  Two independent reasons.
//  (i) At the central element Z = S~^2 = (-I, i) we always have f|Z = i^{-1} f at weight 1/2, so the
//      coset sum is well defined only if rho(Z^{-1}) e_0 = e(1/4) e_0.  WeilRepresentationST gives
//      rho(S)^2 = e(1/4) P (pinned by tests/WeilRepresentation.m) and hence e(-1/4) e_0 -- the wrong
//      sign; the conjugate representation gives e(1/4) e_0.
// (ii) [GY, Lemma 24] attaches q^{-n/M} to the eta with nrd(eta) in n/M + Z, i.e. the exponents of the
//      eta-component are congruent to -Q(eta) mod 1, i.e. rho(T) e_eta = e(-Q(eta)) e_eta.
//
// References: [GY] Guo-Yang; [Sch] Schofer, J. reine angew. Math. 629 (2009); Yifan Yang,
// arXiv:1503.07971.

// ---------------------------------------------------------------------------------------------
// Coset representatives and ST-words
// ---------------------------------------------------------------------------------------------

intrinsic VVCosetReps(M::RngIntElt) -> SeqEnum
{Matrices in SL_2(Z), one per right coset of Gamma_0(M), indexed by P^1(Z/M): Gamma_0(M) g_1 and
 Gamma_0(M) g_2 agree exactly when their bottom rows agree in P^1(Z/M).}
    require M ge 1 : "M must be positive.";
    seen := {};
    reps := [];
    for c in [0..M-1] do
        for d in [0..M-1] do
            if GCD([c, d, M]) ne 1 then continue; end if;
            orb := {@ <(u*c) mod M, (u*d) mod M> : u in [1..M] | GCD(u, M) eq 1 @};
            key := Min([x : x in orb]);
            if key in seen then continue; end if;
            Include(~seen, key);
            // lift (c : d) mod M to a primitive integer pair, then complete the bottom row
            cc := c; dd := d;
            if cc eq 0 and dd eq 0 then cc := M; dd := 1; end if;
            while GCD(cc, dd) ne 1 do dd +:= M; end while;
            _, a, b := XGCD(dd, -cc);                      // a*dd - b*cc = 1
            Append(~reps, Matrix(Integers(), 2, 2, [a, b, cc, dd]));
        end for;
    end for;
    return reps;
end intrinsic;

intrinsic VVSTWord(g::AlgMatElt) -> SeqEnum
{Write g in SL_2(Z) as a word in S = [0,-1;1,0] and T = [1,1;0,1].  Tokens are uniform 2-tuples
 <"S", 0> and <"T", k> (meaning T^k); their product, left to right, is g.  Uniformity matters --
 Magma will not put <"S"> and <"T",k> in one sequence.}
    require Nrows(g) eq 2 and Determinant(g) eq 1 : "Need an element of SL_2(Z).";
    Sm := Matrix(Integers(), 2, 2, [0,-1,1,0]);
    h := g;
    word := [];
    while h[2][1] ne 0 do
        // g = T^k * S * (S^-1 T^-k g), with k chosen so the lower-left entry strictly decreases
        k := Round(h[1][1]/h[2][1]);
        Append(~word, <"T", k>);
        Append(~word, <"S", 0>);
        h := Sm^(-1) * Matrix(Integers(), 2, 2, [1,-k,0,1]) * h;
    end while;
    if h[1][1] eq 1 then
        Append(~word, <"T", h[1][2]>);
    else
        Append(~word, <"S", 0>); Append(~word, <"S", 0>);          // -I = S^2
        Append(~word, <"T", -h[1][2]>);
    end if;
    return word;
end intrinsic;

intrinsic VVWordMatrix(word::SeqEnum) -> AlgMatElt
{The SL_2(Z) matrix of an ST-word, multiplied left to right.}
    Sm := Matrix(Integers(), 2, 2, [0,-1,1,0]);
    h := IdentityMatrix(Integers(), 2);
    for t in word do
        if t[1] eq "S" then h := h*Sm; else h := h*Matrix(Integers(),2,2,[1,t[2],0,1]); end if;
    end for;
    return h;
end intrinsic;

// ---------------------------------------------------------------------------------------------
// Numeric evaluation of eta quotients and of the weight-1/2 metaplectic slash
// ---------------------------------------------------------------------------------------------

intrinsic VVEtaEval(f::EtaQuot, z::FldComElt) -> FldComElt
{f(z), evaluated numerically as sum_r c_r prod_d eta(d z)^r_d.}
    R := Parent(f);
    require Im(z) gt 0 : "z must lie in the upper half plane.";
    etas := [DedekindEta(d*z) : d in R`ds];
    tot := Parent(z)!0;
    for r in Exponents(f) do
        term := Parent(z)!(f`coeffs[r]);
        for i->d in R`ds do if r[i] ne 0 then term *:= etas[i]^r[i]; end if; end for;
        tot +:= term;
    end for;
    return tot;
end intrinsic;

intrinsic VVSlashEval(f::EtaQuot, word::SeqEnum, tau::FldComElt) -> FldComElt
{(f | w)(tau), the weight-1/2 metaplectic slash along the ST-word w.  For w = g_1 g_2 ... g_k one has
 f|w = (...((f|g_1)|g_2)...)|g_k, so evaluating at tau walks the word BACKWARDS, accumulating the
 automorphy factors phi_S(z) = sqrt(z) (principal branch) and phi_T = 1 along the way.}
    z := tau;
    factor := Parent(tau)!1;
    for i := #word to 1 by -1 do
        if word[i][1] eq "S" then
            factor /:= Sqrt(z);
            z := -1/z;
        else
            z := z + word[i][2];
        end if;
    end for;
    return factor * VVEtaEval(f, z);
end intrinsic;

// ---------------------------------------------------------------------------------------------
// The Weil representation over C
// ---------------------------------------------------------------------------------------------

intrinsic WeilRepresentationComplex(Ld::QuaternionLatticeData, CC::FldCom : Dual := true)
      -> AlgMatElt, SeqEnum, SeqEnum, RngIntElt
{rho(S) as a complex matrix, rho(T) as its diagonal, the ordered discriminant-group elements, and the
 index of the trivial coset.  Mirrors WeilRepresentationST, but over ComplexField -- where 1/sqrt|G|
 is unambiguously the positive real root, so the cyclotomic Gauss-sum care of positive_real_sqrt is
 unnecessary.  Dual := true (the default, and what F_f transforms under) conjugates:
     rho*(T) e_eta = e(-nm(eta)) e_eta,   rho*(S)_i,j = e(-1/8)/sqrt|G| * e(+<v_i,v_j>).}
    Q := ChangeRing(Ld`Q, Rationals()); dn := Ld`denom;
    elts := [g : g in Ld`disc_grp]; n := #elts;
    vs := [ChangeRing(g@@Ld`to_disc, Rationals()) : g in elts];
    i0 := rep{i : i in [1..n] | IsZero(elts[i])};
    require &and[&and[IsIntegral(x) : x in Eltseq(v)] : v in vs] :
        "Expected integral lifts of the discriminant group.";

    ii := CC.1; twopii := 2*Pi(CC)*ii;
    sgn := Dual select -1 else 1;
    ee := func<a | Exp(twopii*CC!(a - Floor(a)))>;

    // Gi[i][j] = (v_i Q, v_j) is an integer, and <v_i,v_j> = Gi[i][j]/dn^2, so every exponent lies in
    // (1/L)Z with L = 2 dn^2.  Building the index matrix as ONE compiled mod-L reduction, and then
    // looking roots of unity up in a table, is what keeps |G| = 1800 (3.24e6 entries) tractable.
    Vi := Matrix(Integers(), n, 3, &cat[[Integers()!x : x in Eltseq(v)] : v in vs]);
    Gi := Vi * ChangeRing(Q, Integers()) * Transpose(Vi);
    L := 2*dn^2;
    tab := [ Exp(twopii*CC!(k/L)) : k in [0..L-1] ];
    Gidx := ChangeRing((-sgn*2)*Gi, Integers(L));
    Tidx := ChangeRing(Matrix(Integers(), n, 1, [sgn*Gi[i][i] : i in [1..n]]), Integers(L));

    Tdiag := [ tab[(Integers()!Tidx[i][1]) + 1] : i in [1..n] ];
    c := ee(sgn/8) / Sqrt(CC!n);
    S := ZeroMatrix(CC, n, n);
    for i in [1..n] do
        for j in [1..n] do
            S[i][j] := c * tab[(Integers()!Gidx[i][j]) + 1];
        end for;
    end for;
    return S, Tdiag, elts, i0;
end intrinsic;

intrinsic VVRhoInvE0(S::AlgMatElt, Tdiag::SeqEnum, word::SeqEnum, i0::RngIntElt) -> ModMatFldElt
{rho(gamma^-1) applied to e_0, for gamma given by an ST-word.  Since rho(g_1 ... g_k) is the product
 of the rho(g_i), the inverse applied to e_0 applies rho(g_1)^-1 first.}
    CC := BaseRing(S);
    n := Nrows(S);
    v := ZeroMatrix(CC, n, 1); v[i0][1] := 1;
    for t in word do
        if t[1] eq "S" then
            // rho is unitary and S is symmetric, so S^{-1} = conj(S) and S^{-1} v = conj(S conj(v)).
            // Using this instead of storing a second |G| x |G| matrix halves the memory.
            for i in [1..n] do v[i][1] := ComplexConjugate(v[i][1]); end for;
            v := S * v;
            for i in [1..n] do v[i][1] := ComplexConjugate(v[i][1]); end for;
        else
            for i in [1..n] do v[i][1] *:= Tdiag[i]^(-t[2]); end for;
        end if;
    end for;
    return v;
end intrinsic;

// ---------------------------------------------------------------------------------------------
// The constant terms of F_f
// ---------------------------------------------------------------------------------------------

intrinsic VVConstantTerms(fs::SeqEnum[EtaQuot], Ld::QuaternionLatticeData, M::RngIntElt :
                          Prec := 200, NumSamples := 192, Height := 1, PosDepth := 0)
      -> SeqEnum, SeqEnum, SeqEnum, SeqEnum
{For each f in fs, the constant terms c_eta(0) of F_f at the ISOTROPIC cosets eta.  Returns, per form,
 that sequence of constant terms; then the isotropic cosets themselves (in the matching order, with
 the trivial coset first); then, per form, the maximum deviation of F_f's numerically computed
 PRINCIPAL part from [GY, Lemma 24] in the pipeline's own normalisation.  That last number is the
 correctness gate: it must be small, or the constant terms mean nothing.

 The transcendental work (the coset rho-vectors and the eta values at the pulled-back points) does not
 depend on the form, so passing all the forms of a base at once is far cheaper than one call each.

 With PosDepth = P > 0, a fourth return value gives, per form and per ISOTROPIC coset, the
 POSITIVE coefficients c_eta(j/M) for j = 1..P (so integer index n corresponds to j = n*M).  These
 are the coefficients no scalar q-expansion can supply at a nonzero coset -- the input of the
 weight-3/2 shadow calibration.  ALIASING CAVEAT: the extraction aliases c(nn - K) and c(nn + K),
 so the error at a positive exponent nn grows like exp(4 pi sqrt(p(nn + K)) - 2 pi K + 2 pi nn)
 relative to the constant term; raise NumSamples accordingly and treat the returned gate as
 covering the PRINCIPAL part only.

 CHOOSING THE PARAMETERS -- getting these wrong looks exactly like a wrong answer:
  * Height must be >= 1.  Then Im(gamma tau) <= Im(tau) for every coset, so no coset is pushed up
    into the pole at oo.  Below 1 some coset lands high, individual terms reach 1e65 while their sum
    is O(1), and the cancellation destroys the result.
  * NumSamples must beat the POSITIVE coefficients.  The Fourier extraction aliases c(n - K) and
    c(n + K); for a weakly holomorphic form of pole order p the positive coefficients grow like
    exp(4 pi sqrt(p m)), so the error scales as exp(4 pi sqrt(pK) - 2 pi K Height).  On X0^15(2)'s
    pole-order-30 forms: K = 48 -> 5.8e73, K = 64 -> 1.2e62, K = 96 -> 2.7e28, K = 192 -> 8.7e-103.
    Raising Prec without raising NumSamples does nothing.
  * Prec must cover the dynamic range e^2 pi p Height (1e82 at p = 30) PLUS the cancellation among
    the eta-quotient monomials (about 66 digits on X0^15(2), whose exponent vectors reach +-21).}
    require Height ge 1 : "Height must be at least 1 (see the intrinsic's documentation).";
    CC := ComplexField(Prec);
    ii := CC.1; twopii := 2*Pi(CC)*ii;

    // The S-action goes through the Fourier factorisation (VVWeilFFT): the dense |G| x |G| matrix
    // is what confined this oracle to D*N <= 42, and it is unnecessary.
    fftdata := VVWeilFFT(Ld, CC : Dual := true);
    elts := fftdata[7]; i0 := fftdata[8];
    n := #elts;
    Q := ChangeRing(Ld`Q, Rationals()); dn := Ld`denom;
    vs := [ChangeRing(g@@Ld`to_disc, Rationals()) : g in elts];
    nmv := [ (vs[i]*Q, vs[i])/(2*dn^2) : i in [1..n] ];
    // bucket index j of eta, as in SchoferFormula0 / m0_multiplier; -1 when M*nm is not integral
    res := [ IsIntegral(M*nmv[i]) select (Integers()!(M*nmv[i]) mod M) else -1 : i in [1..n] ];
    iso := [i0] cat [i : i in [1..n] | IsIntegral(nmv[i]) and i ne i0];

    reps := VVCosetReps(M);
    words := [VVSTWord(g) : g in reps];
    U := ZeroMatrix(CC, #reps, n);
    for k in [1..#reps] do
        v := VVRhoInvE0FFT(fftdata, words[k]);
        for i in [1..n] do U[k][i] := v[i]; end for;
    end for;

    // form-independent: the pulled-back point and the eta values there, per (sample, coset)
    KS := NumSamples;
    taus := [CC!(t/KS) + CC!Height*ii : t in [0..KS-1]];
    ds := Parent(fs[1])`ds;
    FACT := ZeroMatrix(CC, KS, #reps);
    ETAB := [[[CC | 0 : d in ds] : k in [1..#reps]] : t in [1..KS]];
    for t in [1..KS] do
        for k in [1..#reps] do
            z := taus[t]; fac := CC!1; w := words[k];
            for r := #w to 1 by -1 do
                if w[r][1] eq "S" then fac /:= Sqrt(z); z := -1/z; else z := z + w[r][2]; end if;
            end for;
            FACT[t][k] := fac;
            ETAB[t][k] := [DedekindEta(d*z) : d in ds];
        end for;
    end for;

    consts := []; errs := []; poscs := [];
    for f in fs do
        require Parent(f)`M eq M : "All forms must live in the eta-quotient ring of level M.";
        foo := qExpansionAtoo(f, 80); f0 := qExpansionAt0(f, 80);
        A := ZeroMatrix(CC, KS, #reps);
        for t in [1..KS] do
            for k in [1..#reps] do
                tot := CC!0;
                for r in Exponents(f) do
                    term := CC!(f`coeffs[r]);
                    for i->d in ds do
                        if r[i] ne 0 then term *:= ETAB[t][k][i]^r[i]; end if;
                    end for;
                    tot +:= term;
                end for;
                A[t][k] := FACT[t][k]*tot;
            end for;
        end for;
        Fv := A*U;                                        // Fv[t][i] = F_{eta_i}(tau_t)
        // Extract every needed Fourier coefficient in ONE matrix product.  Done naively this is
        // |G| x K transcendental evaluations per form (4.6 million at |G| = 24200), which becomes the
        // bottleneck once the S-action is no longer dense; batched it is a single compiled multiply.
        v0 := Minimum(0, Valuation(f0));
        exps := Sort(Setseq( {Rationals() | 0}
                    join {Rationals() | -(Rationals()!j)/M : j in [1..-v0]}
                    join {Rationals() | nn : nn in [Valuation(foo)..-1]}
                    join {Rationals() | (Rationals()!j)/M : j in [1..PosDepth]} ));
        pos := AssociativeArray();
        for a->e in exps do pos[e] := a; end for;
        // NB: build this with an explicit loop.  Magma's multi-index comprehension varies the FIRST
        // index fastest, so [f(a,t) : a in .., t in ..] fills a matrix COLUMN-major -- which silently
        // transposes any non-symmetric matrix built that way.
        E := ZeroMatrix(CC, #exps, KS);
        for a in [1..#exps] do
            for t in [1..KS] do E[a][t] := Exp(-twopii*CC!exps[a]*taus[t]); end for;
        end for;
        Coefs := (E*Fv)/KS;                               // Coefs[a][i] = c_{eta_i}(exps[a])
        coef := func<i, nn | Coefs[pos[Rationals()!nn]][i]>;

        // the gate: F_f's principal part against [GY, Lemma 24] as the pipeline normalises it
        err := 0.0;
        for i in [1..n] do
            for j in [1..-v0] do
                if (j mod M) ne res[i] then continue; end if;
                nn := -(Rationals()!j)/M;
                pred := CC!Coefficient(f0, -j);
                if i eq i0 and IsIntegral(nn) then pred +:= Coefficient(foo, Integers()!nn); end if;
                err := Maximum(err, Abs(coef(i, nn) - pred));
            end for;
        end for;
        for nn in [Valuation(foo)..-1] do
            if (-M*nn) le -v0 then continue; end if;
            err := Maximum(err, Abs(coef(i0, nn) - CC!Coefficient(foo, nn)));
        end for;

        Append(~consts, [coef(i, 0) : i in iso]);
        Append(~errs, err);
        Append(~poscs, [ [coef(i, (Rationals()!j)/M) : j in [1..PosDepth]] : i in iso ]);
    end for;
    return consts, [elts[i] : i in iso], errs, poscs;
end intrinsic;

intrinsic M0MultiplierNumeric(fs::SeqEnum[EtaQuot], Ld::QuaternionLatticeData, D::RngIntElt,
                              N::RngIntElt : Prec := 200, NumSamples := 192, Height := 1)
      -> SeqEnum, SeqEnum
{The m = 0 multipliers of Schofer's formula for the forms fs, computed from the true vector-valued
 F_f, together with the per-form principal-part check errors.  The multiplier is

     (1/2) * c_eta(0)   for any NONZERO ISOTROPIC eta

 (they all carry the same value, and c_0(0) = 0).  Verified against the independently measured
 ground truth on all 9 forms of X0^15(2) and all 5 measured forms of X0^6(5) -- see
 tests/VectorValuedForm.m.  The 1/2 is empirical; it has not been derived.

 This is an ORACLE, not a production route: it costs minutes per base and cannot reach the larger
 discriminant groups.  SchoferFormula.m still uses m0_multiplier; this is what that must reproduce.}
    require N gt 1 : "There are no nonzero isotropic cosets when N = 1, so no m = 0 multiplier.";
    M := IsOdd(D*N) select 4*D*N else 2*D*N;
    consts, isoelts, errs := VVConstantTerms(fs, Ld, M : Prec := Prec, NumSamples := NumSamples,
                                             Height := Height);
    // entry 1 is the trivial coset (whose constant term is 0); any other will do
    return [ c[2]/2 : c in consts ], errs;
end intrinsic;

// ---------------------------------------------------------------------------------------------
// The S-action as a finite Fourier transform on the discriminant group
// ---------------------------------------------------------------------------------------------
// rho(S)_{i,j} = c * e(-sgn <v_i, v_j>) and the pairing is a PERFECT pairing on the finite abelian
// group G = L^v/L, so rho(S) is (up to a relabelling) the Fourier transform on G.  Writing
// G = Z/d_1 x ... x Z/d_k via Moduli, and A_rs = <g_r, g_s> for the generators, put
//         z(y)_r = sum_s (d_r A_rs) y_s   mod d_r      (a bijection of G, by nondegeneracy)
// so that <x, y> = sum_r x_r z(y)_r / d_r.  Then S is a relabelling followed by the standard
// multidimensional DFT, which factors into one small matrix product per axis.
//
// This replaces an |G| x |G| dense matrix by k matrices of sizes d_r x d_r:
//     time   |G|^2            ->  |G| * sum_r d_r
//     memory |G|^2            ->  |G|
// On |G| = 1800 (D*N = 30) that is ~30x fewer operations; on |G| = 24200 (X0^10(11)) it is ~100x,
// and the dense matrix -- 5.9e8 complex entries -- does not have to exist at all.  Without this the
// oracle is confined to D*N <= 42.

intrinsic VVWeilFFT(Ld::QuaternionLatticeData, CC::FldCom : Dual := true) -> List
{Precomputed data for applying rho(S) as a Fourier transform on the discriminant group.  Returns a
 list holding ds, gatherflat, scatterflat, Fmats, c, Tdiag, elts, i0, to be passed to VVApplyS.}
    G := Ld`disc_grp;
    Q := ChangeRing(Ld`Q, Rationals()); dn := Ld`denom;
    mods := Moduli(G);
    keep := [r : r in [1..#mods] | mods[r] gt 1];
    ds := [mods[r] : r in keep];
    k := #ds;
    n := &*ds;
    require n eq #G : "Moduli do not account for the whole discriminant group.";

    elts := [g : g in G];
    require #elts eq n : "Element enumeration disagrees with the group order.";
    i0 := rep{i : i in [1..n] | IsZero(elts[i])};

    // pairing of the generators, A[r][s] = <g_r, g_s> in Q/Z
    gens := [G.(keep[r]) : r in [1..k]];
    wg := [ChangeRing(g@@Ld`to_disc, Rationals()) : g in gens];
    A := [[ (wg[r]*Q, wg[s])/dn^2 : s in [1..k] ] : r in [1..k]];

    strides := [1 : r in [1..k]];
    for r := k-1 to 1 by -1 do strides[r] := strides[r+1]*ds[r+1]; end for;
    flat := func<c | 1 + &+[ (c[r] mod ds[r]) * strides[r] : r in [1..k] ]>;

    coords := [ [ (Eltseq(e)[keep[r]]) mod ds[r] : r in [1..k] ] : e in elts ];
    gatherflat := [ flat([ &+[ Integers()!(ds[r]*A[r][s]) * coords[i][s] : s in [1..k] ]
                           : r in [1..k] ]) : i in [1..n] ];
    scatterflat := [ flat(coords[i]) : i in [1..n] ];
    require #Set(gatherflat) eq n : "The pairing is degenerate: z(y) is not a bijection.";

    ii := CC.1; twopii := 2*Pi(CC)*ii;
    sgn := Dual select -1 else 1;
    Fmats := [* *];
    for r in [1..k] do
        d := ds[r];
        zt := [ Exp(twopii*CC!(-sgn*t/d)) : t in [0..d-1] ];
        Append(~Fmats, Matrix(CC, d, d, [ zt[((x*z) mod d) + 1] : x in [0..d-1], z in [0..d-1] ]));
    end for;
    c := Exp(twopii*CC!(sgn/8)) / Sqrt(CC!n);
    nm := [ (ChangeRing(elts[i]@@Ld`to_disc, Rationals())*Q,
             ChangeRing(elts[i]@@Ld`to_disc, Rationals()))/(2*dn^2) : i in [1..n] ];
    Tdiag := [ Exp(twopii*CC!(sgn*nm[i] - Floor(sgn*nm[i]))) : i in [1..n] ];
    return [* ds, gatherflat, scatterflat, Fmats, c, Tdiag, elts, i0 *];
end intrinsic;

intrinsic VVApplyS(data::List, v::SeqEnum) -> SeqEnum
{Apply rho(S) to a vector indexed in the order of the discriminant-group element list, using the
 Fourier factorisation of VVWeilFFT.}
    ds := data[1]; gatherflat := data[2]; scatterflat := data[3];
    Fmats := data[4]; c := data[5];
    k := #ds; n := &*ds;
    CC := Parent(c);
    w := [CC | 0 : t in [1..n]];
    for i in [1..n] do w[gatherflat[i]] := v[i]; end for;
    // multidimensional DFT: one matrix product per axis, rotating the axes by a transpose so that
    // after k rounds the flat ordering is back where it started.
    W := Matrix(CC, ds[1], n div ds[1], w);
    for r in [1..k] do
        W := Fmats[r] * W;
        nxt := (r eq k) select ds[1] else ds[r+1];
        W := Matrix(CC, nxt, n div nxt, Eltseq(Transpose(W)));
    end for;
    out := Eltseq(W);
    return [ c*out[scatterflat[i]] : i in [1..n] ];
end intrinsic;

intrinsic VVRhoInvE0FFT(data::List, word::SeqEnum) -> SeqEnum
{rho(gamma inverse) applied to e_0 via the Fourier factorisation.  Same convention as VVRhoInvE0:
 the inverse of rho(g_1) is applied first, and S inverse = conj(S) is used (rho is unitary and S is
 symmetric).}
    Tdiag := data[6]; i0 := data[8];
    CC := Parent(data[5]);
    n := #Tdiag;
    v := [CC | 0 : i in [1..n]];  v[i0] := 1;
    for t in word do
        if t[1] eq "S" then
            v := [ComplexConjugate(x) : x in v];
            v := VVApplyS(data, v);
            v := [ComplexConjugate(x) : x in v];
        else
            for i in [1..n] do v[i] *:= Tdiag[i]^(-t[2]); end for;
        end if;
    end for;
    return v;
end intrinsic;
