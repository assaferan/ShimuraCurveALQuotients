// ============================================================================
// ModelVerification.m
//
// Independent sanity checks on the subhyperelliptic cover models stored in
// data/models/models_D_N.m (as produced by AllEquationsAboveCovers / genmodels).
//
// The checks deliberately AVOID the Borcherds/Schofer CM machinery that built the
// models, so they are genuine cross-checks rather than restatements of the same
// computation:
//
//   [1] genus self-consistency : Genus(y^2 + h*y = f) equals the genus recorded
//       next to it in the model file.
//   [2] genus vs theory        : that genus equals X`g for the Shimura quotient
//       X_0(D,N)/W, which comes from the Shimura-curve genus formula.
//   [3] Weil divisibility      : for cover keys W1 subset W2 there is a quotient
//       map X/W1 -> X/W2, so Jac(X/W2) is isogenous to an abelian subvariety of
//       Jac(X/W1); hence over every prime p of good reduction the L-polynomial of
//       X/W2 DIVIDES that of X/W1.  A wrong model generically breaks this.
//
// USAGE (the model file must be `load`-ed first -- Magma's `load` needs a literal):
//   AttachSpec("ShimuraQuotients.spec");
//   load "data/models/models_22_7.m";
//   nchk, nfail := VerifyModelSet(models, 22, 7);
// ============================================================================

function model_lpoly(C, p)
    // numerator of the zeta function of C/F_p; false if p is bad for this model
    try
        Cp := ChangeRing(C, GF(p));
        if not IsNonsingular(Cp) then return false, _; end if;
        return true, Numerator(ZetaFunction(Cp));
    catch e
        return false, _;
    end try;
end function;

// The expected point count #(X_0(D,N)/W)(F_p) comes from ComputePointsViaTrace, which evaluates
// the Eichler-Selberg trace formula on the W-fixed part of the D-new space (TraceDNewALFixed) --
// no modular-symbol space is ever built.  This is the codebase's own trace-formula point count
// (the same routine behind the Weil/automorphism filters), and it is completely independent of
// the Borcherds/Schofer CM machinery that produced the models, so it is a genuine cross-check.

intrinsic VerifyModelSet(models::Assoc, D::RngIntElt, N::RngIntElt : NPrimes := 4, Verbose := true, CheckZeta := true) -> RngIntElt, RngIntElt
{Run independent checks on a model set for X_0(D,N)*. Returns (#checks, #failures).
 CheckZeta runs the trace-formula point-count check.}
    nchk := 0; nfail := 0;
    curves := GetHyperellipticCandidates();

    Cs := AssociativeArray();     // W-set -> model curve
    Xs := AssociativeArray();     // W-set -> the ShimuraQuot X_0(D,N)/W (for the trace-formula check)
    for k in Keys(models) do
        ents := models[k];
        if #ents eq 0 then continue; end if;
        for e in ents do
            if Type(e[2]) eq MonStgElt then continue; end if;   // "CRV": non-hyperelliptic entry
            g := e[1]; f := e[2]; h := e[3];
            C := HyperellipticCurve(f, h);
            gc := Genus(C);
            // [1] genus self-consistency
            nchk +:= 1;
            if gc ne g then
                if Verbose then printf "  [1] FAIL W=%o: recorded genus %o, Genus(model) = %o\n", k, g, gc; end if;
                nfail +:= 1;
            end if;
            // [2] genus vs the Shimura-curve genus formula
            Wset := Set(k);
            if exists(X){Y : Y in curves | Y`D eq D and Y`N eq N and Y`W eq Wset} then
                nchk +:= 1;
                if X`g ne gc then
                    if Verbose then printf "  [2] FAIL W=%o: theory genus %o, model genus %o\n", k, X`g, gc; end if;
                    nfail +:= 1;
                end if;
                Xs[Wset] := X;
            end if;
            Cs[Wset] := C;
        end for;
    end for;
    if Verbose then printf "  model curves: %o\n", #Keys(Cs); end if;

    // [3] Weil-polynomial divisibility across nested cover keys
    Wsets := [W : W in Keys(Cs)];
    pairs := [<W1,W2> : W1 in Wsets, W2 in Wsets | (W1 ne W2) and (W1 subset W2)];
    if #pairs eq 0 then
        // No nested cover pairs (sparse model set): [3] cannot apply -- but [4] still can, and
        // these are exactly the sets that most need it, so fall through rather than returning.
        if Verbose then printf "  [3] no nested cover pairs (skipping [3])\n"; end if;
    end if;
    // primes of good reduction for every model in the set
    ps := []; p := 3;
    while (#pairs gt 0) and (#ps lt NPrimes) and (p lt 200) do
        if (D*N mod p ne 0) then
            ok := true;
            for W in Wsets do
                b := model_lpoly(Cs[W], p);
                if not b then ok := false; break; end if;
            end for;
            if ok then Append(~ps, p); end if;
        end if;
        p := NextPrime(p);
    end while;
    if (#pairs gt 0) and Verbose then printf "  [3] %o nested pair(s); good primes %o\n", #pairs, ps; end if;
    for pr in pairs do
        W1 := pr[1]; W2 := pr[2];            // quotient map X/W1 -> X/W2
        for p in ps do
            b1, L1 := model_lpoly(Cs[W1], p);
            b2, L2 := model_lpoly(Cs[W2], p);
            if not (b1 and b2) then continue; end if;
            nchk +:= 1;
            if not IsDivisibleBy(L1, L2) then
                if Verbose then
                    printf "  [3] FAIL p=%o: L(X/%o) does not divide L(X/%o)\n", p, Sort(SetToSequence(W2)), Sort(SetToSequence(W1));
                end if;
                nfail +:= 1;
            end if;
        end for;
    end for;

    // [4] trace-formula point counts
    if CheckZeta then
        for W in Wsets do
            if not IsDefined(Xs, W) then continue; end if;   // no ShimuraQuot found for this W
            X := Xs[W]; C := Cs[W];
            for p in [3,5,7,11,13,17,19,23] do
                if (D*N mod p eq 0) then continue; end if;
                good := true; actual := 0;
                try
                    Cp := ChangeRing(C, GF(p));
                    if IsNonsingular(Cp) then actual := #Points(Cp); else good := false; end if;
                catch err
                    good := false;
                end try;
                if not good then continue; end if;
                expected := ComputePointsViaTrace(X, p, 1);
                nchk +:= 1;
                if expected ne actual then
                    if Verbose then
                        printf "  [4] FAIL W=%o p=%o: trace formula predicts #X(F_p) = %o, model gives %o\n",
                            Sort(SetToSequence(W)), p, expected, actual;
                    end if;
                    nfail +:= 1;
                end if;
            end for;
        end for;
    end if;
    return nchk, nfail;
end intrinsic;
