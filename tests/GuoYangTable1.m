// EXTERNAL ORACLE: Guo-Yang, "Equations of hyperelliptic Shimura curves",
// Compositio Math. 153 (2017) 1-40, TABLE 1 (p.20) -- "List of hyperelliptic Shimura curves and
// their hyperelliptic involutions", itself quoting Ogg [Ogg83, Theorems 7 and 8].
//
// ⚠ WHY THIS EXISTS. The repo reads Guo-Yang's EQUATION tables (A.1 for level one, A.2 for level
// greater than one) and nothing else. Table 1 is a THIRD published table, on p.20 rather than in the
// appendix, and it was never transcribed. It publishes two columns per base that we can check:
//   * g(D,N), the genus of X_0^D(N);
//   * w, the HYPERELLIPTIC INVOLUTION -- which Atkin-Lehner involution is the hyperelliptic one.
// Found 2026-09-24 by asking of Guo-Yang the question that had just found a whole missing table in
// Gonzalez-Rotger: what else does this paper publish that we do not read?
//
// ⚠ THE SECOND COLUMN IS THE INTERESTING ONE, because w is FREQUENTLY NOT w_{D*N}:
//     57_1 -> w_19    58_1 -> w_29    82_1 -> w_41    93_1 -> w_31
//     10_19 -> w_38   21_2 -> w_7     6_17 -> w_34    10_13 -> w_65
// so it is a real datum about the curve, not a restatement of the level.
//
// ⇒ SHARPNESS, MEASURED BEFORE THIS WAS WRITTEN rather than hoped for. At EVERY ONE of the 43 bases
// the published w is the UNIQUE m | D*N with genus(X_0^D(N)/w_m) = 0. So the check does not merely
// confirm that the published involution is consistent -- it IDENTIFIES it, and this file asserts
// that uniqueness rather than the weaker "the published one qualifies". A test that only checked
// "genus 0" would still pass if the pipeline's genus formula admitted several candidates, which is
// exactly the kind of weak check this repo keeps finding.
//
// COST: the genus formula only -- no Borcherds run, no model read for the genus half. Milliseconds.
//
// ⚠ WHAT THIS DOES NOT DO. It compares our GENUS FORMULA against the published genus, so for the
// first column it is an implementation check against Ogg's list rather than an independent
// derivation. The second column is stronger: which involution is hyperelliptic is a structural fact
// about the curve, and the uniqueness assertion pins it.
SetVerbose("ShimuraQuotients", 0);

printf "Checking Guo-Yang Table 1 (genus and hyperelliptic involution)...";

// <D, N, published g(D,N), published hyperelliptic involution w_m>  -- Table 1, p.20, read
// left-column-then-right-column. 43 rows, which is Ogg's full list of hyperelliptic X_0^D(N), D > 1.
GYT1 := [ <26,1,2,26>,  <35,1,3,35>,  <38,1,2,38>,  <39,1,3,39>,  <51,1,3,51>,  <55,1,3,55>,
          <57,1,3,19>,  <58,1,2,29>,  <62,1,3,62>,  <69,1,3,69>,  <74,1,4,74>,  <82,1,3,41>,
          <86,1,4,86>,  <87,1,5,87>,  <93,1,5,31>,  <94,1,3,94>,  <95,1,7,95>,  <111,1,7,111>,
          <119,1,9,119>, <134,1,6,134>, <146,1,7,146>, <159,1,9,159>, <194,1,9,194>, <206,1,9,206>,
          <6,11,3,66>,  <6,17,3,34>,  <6,19,3,114>, <6,29,5,174>, <6,31,5,186>, <6,37,5,222>,
          <10,11,5,110>, <10,13,3,65>, <10,19,5,38>, <10,23,9,230>,
          <14,3,3,14>,  <14,5,3,14>,  <15,2,3,15>,  <15,4,5,15>,  <21,2,3,7>,
          <22,3,3,66>,  <22,5,5,110>, <26,3,5,26>,  <39,2,7,39> ];
error if #GYT1 ne 43,
    Sprintf("GuoYangTable1: %o rows, not the 43 of Ogg's list -- fix the transcription before "
            * "trusting any verdict below", #GYT1);

n_genus := 0; n_invol := 0; gbad := []; wbad := []; sbad := [];
for t in GYT1 do
    D := t[1]; N := t[2]; gpub := t[3]; mpub := t[4];

    // --- column 1: the genus ---
    gours := GenusShimuraCurveQuotient(D, N, {Integers()|1});
    n_genus +:= 1;
    if gours ne gpub then
        Append(~gbad, Sprintf("%o_%o: ours g = %o, Guo-Yang Table 1 publishes %o", D, N, gours, gpub));
    end if;

    // --- column 2: the hyperelliptic involution, asserted by UNIQUENESS ---
    // w_m is the hyperelliptic involution exactly when X/w_m has genus 0 (the quotient of a
    // hyperelliptic curve by its hyperelliptic involution is P^1). Collect every candidate.
    cands := [m : m in Divisors(D*N) | m gt 1 and GCD(m, (D*N) div m) eq 1
                  and GenusShimuraCurveQuotient(D, N, {Integers()|1, m}) eq 0];
    n_invol +:= 1;
    if mpub notin cands then
        Append(~wbad, Sprintf("%o_%o: Guo-Yang publish w_%o as the hyperelliptic involution, but "
                              * "X/w_%o has genus %o, not 0 (genus-0 quotients here: %o)",
                              D, N, mpub, mpub,
                              GenusShimuraCurveQuotient(D, N, {Integers()|1, mpub}), cands));
    elif #cands ne 1 then
        // ⚠ NOT a disagreement with the paper -- the published w is still among the candidates.
        // It means the check has STOPPED BEING SHARP at this base, so it no longer identifies w.
        // Measured 2026-09-24: #cands = 1 at all 43. If this fires, the genus formula changed.
        Append(~sbad, Sprintf("%o_%o: w_%o is published and does qualify, but %o involutions give a "
                              * "genus-0 quotient (%o), so this base no longer PINS the "
                              * "hyperelliptic involution", D, N, mpub, #cands, cands));
    end if;
end for;

error if not IsEmpty(gbad),
    Sprintf("Guo-Yang Table 1: %o genus disagreement(s) with the published list: %o", #gbad, gbad);
error if not IsEmpty(wbad),
    Sprintf("Guo-Yang Table 1: %o hyperelliptic-involution disagreement(s): %o", #wbad, wbad);
error if not IsEmpty(sbad),
    Sprintf("Guo-Yang Table 1: %o base(s) where the published involution is no longer UNIQUELY "
            * "pinned: %o.\n  This is not a conflict with the paper -- it means the check weakened "
            * "from an identification to a consistency check, which is worth knowing.", #sbad, sbad);
// ⚠ NON-VACUITY: 43 rows, so both counters must reach 43. A transcription that silently lost rows,
// or a loop that stopped executing, would otherwise read as a clean pass.
error if n_genus ne 43 or n_invol ne 43,
    Sprintf("Guo-Yang Table 1: %o genus and %o involution comparison(s), expected 43 of each -- "
            * "something stopped being compared", n_genus, n_invol);

printf " ok (%o genus + %o hyperelliptic-involution check(s) against Guo-Yang Table 1; the published "
       * "involution is the UNIQUE genus-0 quotient at every base)\n", n_genus, n_invol;
