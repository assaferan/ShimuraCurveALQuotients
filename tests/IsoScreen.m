// tests/IsoScreen.m
//
// A FAST, SOUND SCREEN for "is our model the same curve as the published one?"
//
// WHY. Verifying a new base against Guo-Yang means an isomorphism test, and the cost of
// IsIsomorphic depends violently on the PRESENTATION, not just the genus:
//
//     39_2  W={1}  hyperelliptic, genus 7   IsIsomorphic in 0.06 s
//     21_2  W={1}  CRV pair,      genus 3   ~100 s
//     14_3  W={1}  CRV pair,      genus 3   still running after 50+ min
//     26_3  W={1}  CRV pair,      genus 5   still running after 1 h
//
// So a plain IsIsomorphic is unusable as a first look at a complete intersection in weighted
// projective space. This screens in seconds per prime instead, and is decisive on a mismatch.
//
// THE INVARIANT. Count degree-1 PLACES of the function field over GF(p). That is the smooth
// projective point count, and it IS an isomorphism invariant.
//
// ⚠⚠ DO NOT SCREEN WITH AFFINE POINT COUNTS. They depend on the affine chart, so isomorphic
// curves can disagree, and a mismatch proves nothing. This is not hypothetical: on 2026-09-05 an
// affine screen of 14_3 reported mismatches at p = 13, 29, 53 (12 vs 16, 44 vs 40, 36 vs 32) and
// would have condemned a model that the place counts then matched at ALL THIRTEEN primes. Every
// one of those three was a chart artifact. Affine counts nearly produced a false NEGATIVE, which
// is worse than a false positive here -- it discards a correct result.
//
// ⚠ Also do not screen with #Points on a curve in a WeightedProjectiveSpace: it returned nothing
// at every prime tried (silently, so the run looked like agreement until a guard was added).
//
// HOW TO READ THE RESULT.
//   * a MISMATCH at one prime disproves isomorphism outright;
//   * agreement over many primes is STRONG EVIDENCE, not proof -- non-isomorphic curves can share
//     point counts. Before committing a model, still get an exact IsIsomorphic (a hyperelliptic
//     quotient is often tractable when the full curve is not: 26_3's genus-2 [1,78] cover settled
//     in seconds where its genus-5 W={1} did not).
//   * ncmp = 0 means the test was VACUOUS -- report that, never "no mismatches".

// Count degree-1 places of C over GF(p); returns <ok, n>.
function places_over(C_at_p)
    try
        K := AlgorithmicFunctionField(FunctionField(C_at_p));
        return true, #Places(K, 1);
    catch e
        return false, 0;
    end try;
end function;

// build_o / build_g take (F, S, Y, X) and return defining equations in AffineSpace(F,3).
// badp: primes to skip (bad reduction / denominators) -- pass 2*D*N times any extra.
// Returns <ncompared, nmismatch>, and prints a per-prime table.
function ScreenByPlaces(name, badp, build_o, build_g, primes)
    printf "IsoScreen %o:  p | ours | published | match\n", name;
    ncmp := 0; nmis := 0;
    for p in primes do
        if badp mod p eq 0 then continue; end if;
        F := GF(p);
        A<S,Y,X> := AffineSpace(F, 3);
        ok1 := false; ok2 := false; n1 := 0; n2 := 0;
        try
            ok1, n1 := places_over(Curve(A, build_o(F, S, Y, X)));
            ok2, n2 := places_over(Curve(A, build_g(F, S, Y, X)));
        catch e
            ok1 := false;
        end try;
        if not (ok1 and ok2) then printf "  %o | (reduction failed)\n", p; continue; end if;
        ncmp +:= 1;
        if n1 ne n2 then nmis +:= 1; end if;
        printf "  %o | %o | %o | %o\n", p, n1, n2, n1 eq n2;
    end for;
    // Report the COMPARISON COUNT, not just the mismatch count: a run where every reduction
    // failed has zero mismatches and must not be mistaken for agreement.
    printf "  compared %o prime(s), %o mismatch(es) -> ", ncmp, nmis;
    if ncmp eq 0 then
        printf "NO EVIDENCE (vacuous -- every reduction failed)\n";
    elif nmis gt 0 then
        printf "NOT ISOMORPHIC (a degree-1 place count is an invariant, so one mismatch settles it)\n";
    else
        printf "consistent over %o primes (strong evidence, NOT proof)\n", ncmp;
    end if;
    return ncmp, nmis;
end function;

// Nothing runs at load: this file is a helper, imported as
//     import "tests/IsoScreen.m" : ScreenByPlaces;
// Worked example (14_3 against Guo-Yang's z^2=-9x^2-2, y^2=-7x^4+22x^2+1), 13/13 consistent:
//   ScreenByPlaces("14_3", 2*7*3,
//     func<F,S,Y,X | [ Y^2 + (F!256*S^4 + (F!1792/3)*S^3 + (F!6592/9)*S^2 + (F!1792/3)*S + F!256),
//                      X^2 + (F!3*S^2 + F!2*S + F!3) ]>,
//     func<F,S,Y,X | [ Y^2 - (-F!7*S^4 + F!22*S^2 + 1), X^2 - (-F!9*S^2 - 2) ]>,
//     [5,11,13,17,19,23,29,31,37,41,43,47,53]);
