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

intrinsic M0MultiplierExact(fs::SeqEnum[EtaQuot], Ld::QuaternionLatticeData, D::RngIntElt,
                            N::RngIntElt : Prec := 80) -> SeqEnum
{The m = 0 multipliers of Schofer's formula for the forms fs, evaluated EXACTLY as a finite sum
 over Gamma_0(M) cosets -- no Fourier sampling, no CM table, minutes per base.  The multiplier is

     (1/2) * c_eta(0)   for any NONZERO ISOTROPIC eta

 (all such eta carry the same value, and c_0(0) = 0), and c_eta(0) is assembled coset by coset:
 for each coset word w and eta monomial prod_d eta(d tau)^(r_d), triangularising
 [d 0; 0 1] g = g_d [a_d b_d; 0 e_d] exhibits f|w as a constant multiple of a q-series with
 exact rational exponents and root-of-unity coefficients; the constant per (monomial, word) is
 pinned numerically at one point of the upper half plane and verified at a second (its closed
 form -- the Dedekind-sum eta multiplier system -- is certified in vvdata/weyl-campaign/cusp4.m).
 Then c_eta(0) = sum_w rho(w^-1)e_0[eta] * a0(f|w).

 Validated against the measured ground truth on 21 bases (vvdata/weyl-campaign, branch
 m0-theta-campaign): the 15_2 full panel 9/9 exactly, 21_2, 22_3, and every base of the
 constraints ledger.  Returns one rational per form; raises an error rather than return an
 unverified value (kappa two-point check, isotropic-component agreement, rational snap).}
    require N gt 1 : "There are no nonzero isotropic cosets when N = 1, so no m = 0 multiplier.";
    require #fs gt 0 : "Empty form sequence.";
    R := Parent(fs[1]); ds := R`ds;
    M := IsOdd(D*N) select 4*D*N else 2*D*N;
    require R`M eq M : "The forms live at the wrong level for (D, N).";

    CC := ComplexField(Prec); ii := CC.1; pi := Pi(CC);
    ee := func< z | Exp(2*pi*ii*z) >;

    fftdata := VVWeilFFT(Ld, CC : Dual := true);
    elts := fftdata[7]; i0 := fftdata[8];
    Qr := ChangeRing(Ld`Q, Rationals()); dn := Ld`denom;
    isoidx := [];
    for i in [1..#elts] do
        v := ChangeRing(elts[i]@@Ld`to_disc, Rationals());
        r := (v*Qr, v)/(2*dn^2);
        if r eq Floor(r) then Append(~isoidx, i); end if;
    end for;
    require #isoidx ge 2 : "No nonzero isotropic coset found.";

    reps := VVCosetReps(M);
    words := [ VVSTWord(g) : g in reps ];

    triang := function(g, d)
        g2 := Matrix(Integers(), 2, 2, [d*g[1][1], d*g[1][2], g[2][1], g[2][2]]);
        c1 := g2[1][1]; c2 := g2[2][1];
        h := GCD(c1, c2);
        p1 := c1 div h; p2 := c2 div h;
        gg, u, v := XGCD(p1, p2);
        error if gg ne 1, "triangularisation: row not primitive";
        gd := Matrix(Integers(), 2, 2, [p1, -v, p2, u]);
        sd := gd^(-1) * g2;
        a := sd[1][1]; b := sd[1][2]; e := sd[2][2];
        if a lt 0 then a := -a; b := -b; e := -e; end if;
        if e lt 0 then gd := -gd; sd := -sd; a := sd[1][1]; b := sd[1][2]; e := sd[2][2]; end if;
        return a, b, e;
    end function;

    slashdata := function(word, tau)
        z := tau; factor := CC!1;
        for i := #word to 1 by -1 do
            if word[i][1] eq "S" then factor /:= Sqrt(z); z := -1/z;
            else z := z + word[i][2]; end if;
        end for;
        return factor, z;
    end function;

    tau0 := CC!0.31 + CC!1.31*ii;
    tau1 := CC!(-0.57) + CC!1.73*ii;

    // PER-WORD FALLBACK for the two evaluation points, triggered only when the global tau0/tau1
    // land a word's z = w.tau outside the measured-safe range. eta(d z) is evaluated for d up to
    // M, so Im(z) too small OR too large both blow up its dynamic range in floating point; the
    // measured bands (vvdata/weyl-campaign/tau-precision/tauwindow.m, over all 2232 cosets of
    // M = 1220) show catastrophic loss ONLY at the extremes -- [0,1e-5) and [1e-1,1e2) lose
    // 90-180 digits, everything from [1e-5,1e-1) loses at most ~10, comfortably inside the
    // 1e-15-relative two-point check below (needs ~15 valid digits out of Prec := 80).
    //
    // THE UPPER THRESHOLD SCALES WITH M, the lower one does not -- and using the raw M = 1220
    // absolute band on a smaller base is a real bug (caught by M0PROGRESS=1 on 15_2, M = 60: it
    // flagged the DEFAULT tau0 itself, Im(z) = 1.31, as unsafe, when 15_2 has always worked
    // exactly with that tau0). The mechanism explains why: eta is evaluated at d*z for d up to
    // M, so the dynamic-range risk on the large side is really about M*Im(z), not Im(z) alone --
    // the M = 1220 measurement's unsafe boundary sits at Im(z) ~ 0.1, i.e. M*Im(z) ~ 122. The
    // near-real-axis mechanism on the small side (z and -1/z both close to the real line) is
    // about the word's own S-step structure, not a d-scaling effect, so that threshold is kept
    // as an absolute floor.
    //
    // The correction: since the slash constant is tau-independent (that is exactly what the
    // two-point check verifies), tau is a free choice per word. Take tau to be the PREIMAGE of a
    // fixed safe-band target z under the word's own SL2(Z) matrix (tau := g^-1.ztarg), so
    // slashdata(w, tau) reproduces z = ztarg exactly regardless of how extreme that word's own
    // (c, d) are. This is deliberately NOT applied to every word: an earlier attempt that used
    // per-word preimages everywhere pushed Im(tau) itself to extreme values for many words
    // (sfun's own argument is evaluated at tau directly, not at z), which was slow rather than
    // wrong but never finished on 58_5 in a reasonable time. Confining it to the rare words the
    // default actually mishandles keeps that risk rare too.
    z0targ := CC!0.31 + CC!0.0057*ii;
    z1targ := CC!(-0.57) + CC!0.0083*ii;
    safe_im := func< z | Abs(Im(z)) ge 10^(-5) and M*Abs(Im(z)) lt 100 >;
    preimage := function(g, z)
        a := CC!g[1][1]; b := CC!g[1][2]; c := CC!g[2][1]; d := CC!g[2][2];
        return (d*z - b) / (a - c*z);
    end function;

    monos := {@ @};
    for f in fs do for r in Exponents(f) do Include(~monos, r); end for; end for;

    // Cusp-class shortcut: the per-coset product rho(w^-1)e_0[eta] * a0(f|w) is CONSTANT
    // on each class g = gcd(c, M) for every isotropic component (e_0 is a
    // rho(Gamma_0(lev))-eigenvector whose character is f's eta multiplier, and the
    // isotropic components are T-invariant; verified on 570 class/form/base checks and
    // derived in vvdata/weyl-campaign/thetag-derivation.md).  So evaluate ONE canonical
    // coset per class, VERIFY the constancy on up to two more cosets of the class, and
    // multiply by the true class size.  This replaces #cosets FFT applications by
    // ~3 * #classes -- the difference between minutes and hours at M ~ 400.
    classof := [ GCD(VVWordMatrix(w)[2][1] mod M, M) : w in words ];
    classes := Sort(Setseq(Set(classof)));
    Ng := AssociativeArray();
    canon := AssociativeArray();   // class -> [canonical wi, up to 2 verification wi's]
    for g0 in classes do
        idxs := [ wi : wi in [1..#words] | classof[wi] eq g0 ];
        Ng[g0] := #idxs;
        picks := [ idxs[1] ];
        if #idxs ge 2 then Append(~picks, idxs[1 + (#idxs div 2)]); end if;
        if #idxs ge 3 then Append(~picks, idxs[#idxs]); end if;
        canon[g0] := picks;
    end for;
    selected := Sort(Setseq(&join{ Set(canon[g0]) : g0 in classes }));

    // Diagnostic-only, gated behind M0PROGRESS so it's silent by default: printf is buffered to
    // a file (CLAUDE.md), so a killed/crashed run leaves NOTHING to show where it got to.
    // WriteStderr is not buffered that way.
    progress := GetEnv("M0PROGRESS") ne "";
    if progress then
        WriteStderr(Sprintf("M0MultiplierExact: %o words selected of %o total, %o classes\n",
                             #selected, #words, #classes));
    end if;

    SS := PowerSeriesRing(CC); t := SS.1;
    a0tab := [ [ CC!0 : r in monos ] : w in words ];
    nfallback := 0;
    t_progress := Realtime();
    for wcount->wi in selected do
        w := words[wi];
        g := VVWordMatrix(w);
        tri := [ ];
        for d in ds do
            a, b, e := triang(g, d);
            Append(~tri, <a, b, e>);
        end for;
        W := LCM([ tri[i][3] : i in [1..#ds] ]);
        leads := [ &+[ Integers() | r[i]*tri[i][1]*(W div tri[i][3]) : i in [1..#ds] ]
                   : r in monos ];
        depth := Maximum([ 0 ] cat [ -L : L in leads ]) + 1;
        if progress then
            WriteStderr(Sprintf("  wi=%o W=%o depth=%o\n", wi, W, depth));
        end if;
        units := [ ];
        for i->d in ds do
            a, b, e := Explode(tri[i]);
            step := 24*a*(W div e);
            u := SS!1 + O(t^depth);
            n := 1;
            while n*step lt depth do
                u *:= 1 - ee(CC!(n*b/e))*t^(n*step);
                n +:= 1;
            end while;
            Append(~units, u);
        end for;
        mytau0 := tau0; mytau1 := tau1;
        _, ztest0 := slashdata(w, tau0);
        if not safe_im(ztest0) then
            mytau0 := preimage(g, z0targ); nfallback +:= 1;
            if progress then
                WriteStderr(Sprintf("  fallback wi=%o point0: default Im(z)=%o -> using preimage\n",
                                     wi, RealField(6)!Im(ztest0)));
            end if;
        end if;
        _, ztest1 := slashdata(w, tau1);
        if not safe_im(ztest1) then
            mytau1 := preimage(g, z1targ); nfallback +:= 1;
            if progress then
                WriteStderr(Sprintf("  fallback wi=%o point1: default Im(z)=%o -> using preimage\n",
                                     wi, RealField(6)!Im(ztest1)));
            end if;
        end if;
        fac0, z0 := slashdata(w, mytau0);
        fac1, z1 := slashdata(w, mytau1);
        if progress and (wcount mod 20 eq 0 or Realtime(t_progress) gt 30) then
            WriteStderr(Sprintf("  word %o/%o (wi=%o), %o fallbacks so far\n",
                                 wcount, #selected, wi, nfallback));
            t_progress := Realtime();
        end if;
        for ri->r in monos do
            if progress and ri mod 50 eq 0 then
                WriteStderr(Sprintf("    wi=%o mono %o/%o\n", wi, ri, #monos));
            end if;
            L := leads[ri];
            if L gt 0 then continue; end if;
            produ := &*[ SS | units[i]^(r[i]) : i in [1..#ds] | r[i] ne 0 ];
            c0 := Coefficient(produ, -L);
            if c0 eq 0 then a0tab[wi][ri] := CC!0; continue; end if;
            num0 := fac0 * &*[ CC | DedekindEta(d*z0)^(r[i]) : i->d in ds | r[i] ne 0 ];
            num1 := fac1 * &*[ CC | DedekindEta(d*z1)^(r[i]) : i->d in ds | r[i] ne 0 ];
            sfun := func< tau | ee(tau*L/(24*W)) *
                &*[ CC | ( DedekindEta((tri[i][1]*tau + tri[i][2])/tri[i][3]) *
                           ee(-(tri[i][1]*tau + tri[i][2])/(24*tri[i][3])) )^(r[i])
                    : i in [1..#ds] | r[i] ne 0 ] >;
            k0 := num0 / sfun(mytau0);
            k1 := num1 / sfun(mytau1);
            // The two-point check must SCALE WITH THE CONSTANT, at the same 10^(-15) the other
            // four guards in this intrinsic use (10^(-15)*vscale on the value checks,
            // 10^(-15)*Max(1,|.|) on the contributions).  This was the last absolute tolerance
            // here; the merged m0exact-relative-tolerance work converted the others but not this.
            //
            // MEASURED, two bases, reporting instead of erroring:
            //   X0^58(5)   1446 evaluations   |k| 1 .. 2.5e10   reldiff <= 5.9e-33
            //   X0^34(11)   319 failures      |k| 1e8 .. 1e12   reldiff <= 1.07e-18
            // The absolute 1e-30 failed on LARGE constants whose agreement was unchanged
            // (absdiff = reldiff * |k|), which is why it had to become relative.  But the
            // ACHIEVABLE relative precision is base-dependent: 33 digits at 58_5, 18 at 34_11,
            // at the same Prec := 80 -- longer eta products accumulate more rounding.  So the
            // threshold must be the siblings' 10^(-15), not something tuned to the best base.
            // (A first attempt used 10^(-30) relative, calibrated on 58_5 alone; that passed
            // 58_5 and still blocked 34_11 and 74_5.  Do not re-tighten it on one base.)
            //
            // 10^(-15) still catches what this guard is for: a wrong slash constant differs at
            // O(1), not in the 19th digit.  Maximum(1, .) keeps it exactly as strict as the
            // original absolute test for |k| <= 1.
            kscale := Maximum(Abs(k0), Abs(k1));
            error if Abs(k0 - k1) gt 10^(-15) * Maximum(1, kscale),
                "M0MultiplierExact: slash constant failed its two-point check";
            a0tab[wi][ri] := k0 * c0;
        end for;
    end for;
    if progress then
        WriteStderr(Sprintf("M0MultiplierExact: a0 table done, %o fallback points of %o\n",
                             nfallback, 2*#selected));
    end if;

    rvtab := AssociativeArray();
    for wi in selected do rvtab[wi] := VVRhoInvE0FFT(fftdata, words[wi]); end for;

    mults := [ Rationals() | ];
    for f in fs do
        a0w := AssociativeArray();
        for wi in selected do
            a0w[wi] := &+[ CC | f`coeffs[r] * a0tab[wi][Index(monos, r)]
                           : r in Exponents(f) ];
        end for;
        cvals := [ CC!0 : i in isoidx ];
        for g0 in classes do
            picks := canon[g0];
            contribs := [ [ rvtab[wi][i] * a0w[wi] : i in isoidx ] : wi in picks ];
            // constancy check across the class's sampled cosets, per component.
            // RELATIVE tolerance with a 1e-15 floor: 39_2 (M = 156, the deepest
            // base, fractional-pp family) shows measured class-1 deviation
            // 1.1e-22 at Prec 80 -- honest deep-series roundoff, while its
            // cusp3 dump passes every class-constancy check (cuspclass2.py,
            // 42/42), so constancy itself is exact.  A genuine violation is
            // O(scale), so 1e-15 relative keeps 15 digits of margin.  Max(1, .)
            // keeps an absolute floor at small scales.
            for k := 2 to #picks do
                for j := 1 to #isoidx do
                    dev := Abs(contribs[k][j] - contribs[1][j]);
                    tol := 10^(-15) * Maximum(1, Abs(contribs[1][j]));
                    error if dev gt tol,
                        Sprintf("M0MultiplierExact: class-constancy violated (class %o, dev %o, scale %o)",
                                g0, RealField(6)!dev, RealField(6)!Abs(contribs[1][j]));
                end for;
            end for;
            for j := 1 to #isoidx do
                cvals[j] +:= Ng[g0] * contribs[1][j];
            end for;
        end for;
        // the nonzero isotropic components must agree, be real, and snap to a rational
        // (same relative-tolerance reasoning as the constancy check above)
        vals := [ cvals[j] : j->i in isoidx | i ne i0 ];
        vscale := Maximum(1, Abs(vals[1]));
        for v in vals do
            error if Abs(v - vals[1]) gt 10^(-15) * vscale,
                Sprintf("M0MultiplierExact: isotropic components disagree (dev %o, scale %o)",
                        RealField(6)!Abs(v - vals[1]), RealField(6)!vscale);
        end for;
        error if Abs(Im(vals[1])) gt 10^(-15) * vscale,
            "M0MultiplierExact: constant term not real";
        mult := BestApproximation(Re(vals[1])/2, 10^4);
        error if Abs(CC!mult - Re(vals[1])/2) gt 10^(-15) * vscale,
            "M0MultiplierExact: multiplier does not snap to a rational";
        Append(~mults, mult);
    end for;
    return mults;
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
