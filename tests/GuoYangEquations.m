// Do our committed models actually reproduce Guo-Yang's PUBLISHED equations?
//
// WHY THIS EXISTS. `ModelChecks.m` validates every model file structurally -- genus
// self-consistency, the Shimura-curve genus formula, Weil-polynomial divisibility, trace-formula
// point counts -- but it never compares a model to the literature. So nothing in CI would notice
// a committed model silently disagreeing with the paper it is supposed to reproduce. The
// `X0_D_N.m` tests do compare against Guo-Yang, but they re-run `AllEquationsAboveCovers` (minutes
// to hours each), which is why only some bases have one. This test is the cheap half: it reads the
// STORED model and compares curves directly, so it costs milliseconds and can cover any base.
//
// SOURCE. J.-W. Guo and Y. Yang, "Equations of hyperelliptic Shimura curves" (arXiv:1510.06193),
// the "Equations of level ..." tables.
//
// ⚠ EVERY ENTRY BELOW IS HAND-TRANSCRIBED, DELIBERATELY. Automated extraction from the arXiv
// LaTeX was tried on 2026-09-04 and abandoned; see PLAN.md, COVERAGE. It hit four separate
// silent-corruption bugs (equations wrapping across `\\` into several $...$ groups; a leading
// minus belonging to the first term rather than being an overall factor; $-parity broken by
// segmenting at the label, which is itself inside $...$; newlines inside the math), each of which
// produces a plausible WRONG polynomial rather than an error. The tables are also heterogeneous:
// some rows are a single y^2=f(x), others a PAIR (82_1: y^2=f(s) AND x^2=g(s)), 15_4 is a conic in
// z, and 93_1 mixes two variables. Do not re-automate this.
// ⚠ UPDATE 2026-09-05: that 93_1 anomaly is now SETTLED, not merely suspected. Their `-3t` is a
// typo for `-3s` -- our model's genus-2 quotient is isomorphic to that reading and to none of the
// three other plausible repairs. See the 93_1 block at the end of this file. It is a worked example
// of why hand-transcription is the right call: an extractor would have emitted a two-variable
// polynomial or silently dropped a term, and either way produced a wrong curve without an error.
//
// THE PUBLISHED EQUATION CORRESPONDS TO OUR W={1} KEY (the curve itself). Each case still records
// the key explicitly, but the mapping is uniform as far as checked.
// ⚠ A CORRECTION, recorded because the wrong version briefly stood: an earlier note here claimed
// `X0_26_1.m` compares a degree-3 s-equation against Guo-Yang's degree-6 x-equation, and inferred
// the mapping was not uniform. That was an artifact of grepping only the FIRST HyperellipticCurve
// in the file. X0_26_1.m in fact carries THREE cover entries, and its {1} entry is exactly
// Guo-Yang's -2s^6+19s^4-24s^2-169 with the identity matrix.
//
// HOW THE PAIRED / CONIC PRESENTATIONS ARE HANDLED (see X0_82_1.m, X0_26_1.m, X0_15_1.m,
// X0_146_1.m): build the pair as a Curve in a WeightedProjectiveSpace -- e.g. weights [1,2,1,1]
// for (x,y,z,s), since y has weight 2 -- and pin the isomorphism by giving cover_data an entry for
// EACH cover, the quotients as well as the full curve, each with its own matrix. Establishing the
// quotients first is what makes the full-curve isomorphism findable.
//
// COST: the hyperelliptic comparisons are ~0.03s each; the single PAIRED case (21_2) is a genus-3
// complete intersection and IsIsomorphic there takes ~100s, dominating the file. Still fine for
// CI, but do not assume more paired bases are free.
//
// The comparison is `IsIsomorphic` on the curves, not coefficient equality: a model is only
// defined up to isomorphism, and ours generally differs from the published one by a coordinate
// change (51_1's is x -> 3x, y -> (27/16)y).

// NB: TOP-LEVEL statements, not a procedure. run_tests.m runs tests via `eval`, and an `eval`
// inside a procedure that closes over an outer variable SEGFAULTS Magma 2.29 -- the same trap
// ModelChecks.m documents, and which this file hit on first writing.
printf "Comparing committed models against Guo-Yang's published equations...";
P<x> := PolynomialRing(Rationals());

// <D, N, cover key (Sort(W)), published f with y^2 = f(x)>
gy_cases := [*
    // X_0^51(1):  y^2 = -(x^2+3)(243x^6+235x^4-31x^2+1)
    // ours is the W={1} genus-3 entry; matches under x -> 3x, y -> (27/16)y (over-determined,
    // fixed by the leading and constant coefficients and consistent on all five).
    <51, 1, [Integers()|1], -(x^2+3)*(243*x^6+235*x^4-31*x^2+1)>,

    // X_0^87(1):  y^2 = -(x^6-7x^4+43x^2+27)(243x^6+523x^4+369x^2+81)   [genus 5]
    <87, 1, [Integers()|1], -(x^6-7*x^4+43*x^2+27)*(243*x^6+523*x^4+369*x^2+81)>,

    // X_0^14(5):  y^2 = -23x^8 - 180x^7 - 358x^6 - 168x^5 - 677x^4 + 168x^3 - 358x^2 + 180x - 23
    // ⚠ In the source this equation WRAPS across `\\` into two $...$ groups; the tail
    // (+168x^3 - 358x^2 + 180x - 23) is easy to lose. Transcribed whole.                 [genus 3]
    <14, 5, [Integers()|1],
      -23*x^8 - 180*x^7 - 358*x^6 - 168*x^5 - 677*x^4 + 168*x^3 - 358*x^2 + 180*x - 23>,

    // X_0^55(1):  y^2 = -(x^4-x^3+x^2+x+1)(3x^4+x^3-5x^2-x+3)                            [genus 3]
    <55, 1, [Integers()|1], -(x^4-x^3+x^2+x+1)*(3*x^4+x^3-5*x^2-x+3)>,

    // X_0^15(2):  y^2 = -(x^2+3)(3x^2+4)(x^4-x^2+4)                                      [genus 3]
    <15, 2, [Integers()|1], -(x^2+3)*(3*x^2+4)*(x^4-x^2+4)>,

    // X_0^22(3):  y^2 = -27x^8 - 308x^6 - 2146x^4 - 308x^2 - 27                          [genus 3]
    <22, 3, [Integers()|1], -27*x^8 - 308*x^6 - 2146*x^4 - 308*x^2 - 27>,

    // X_0^39(2):  y^2 = -(x^8+11x^7+52x^6+140x^5+243x^4+280x^3+208x^2+88x+16)
    //                   (7x^4+24x^3+32x^2+24x+16)(x^4+3x^3+8x^2+12x+7)          [degree 16, genus 7]
    // ⚠ In the source this WRAPS across `\\`, the second line beginning `\qquad(7x^4...` -- reading
    // only the first math group gives a degree-8 factor and a plausible WRONG curve.
    // ⚠ models_39_2.m was produced with CMNONCOPRIME=1 and does NOT regenerate by default; see its
    // header. This comparison is the independent oracle that makes that file trustworthy, so it is
    // the one test that must not be skipped when the CM filter is changed.
    <39, 2, [Integers()|1],
      -(x^8+11*x^7+52*x^6+140*x^5+243*x^4+280*x^3+208*x^2+88*x+16)
       *(7*x^4+24*x^3+32*x^2+24*x+16)*(x^4+3*x^3+8*x^2+12*x+7)>
*];

gy_checked := 0;
for c in gy_cases do
    gy_D, gy_N, gy_key, gy_f := Explode(c);
    gy_mf := Sprintf("data/models/models_%o_%o.m", gy_D, gy_N);
    gy_models := eval (Read(gy_mf) cat "\nreturn models;");
    gy_ok, gy_entry := IsDefined(gy_models, gy_key);
    error if not gy_ok,
        Sprintf("X0^%o(%o): model file has no cover key %o", gy_D, gy_N, gy_key);
    gy_g, gy_fo, gy_h := Explode(gy_entry[1]);
    error if Type(gy_fo) eq MonStgElt,
        Sprintf("X0^%o(%o): cover %o is stored as a general CRV, not a hyperelliptic model; "
                * "it needs a different comparison", gy_D, gy_N, gy_key);
    gy_Cours := HyperellipticCurve(gy_fo, gy_h);
    gy_Cgy   := HyperellipticCurve(gy_f);
    error if Genus(gy_Cours) ne Genus(gy_Cgy),
        Sprintf("X0^%o(%o) cover %o: our genus %o vs Guo-Yang's %o -- comparing the wrong object",
                gy_D, gy_N, gy_key, Genus(gy_Cours), Genus(gy_Cgy));
    error if not IsIsomorphic(gy_Cours, gy_Cgy),
        Sprintf("X0^%o(%o) cover %o: our model is NOT isomorphic to Guo-Yang's published curve",
                gy_D, gy_N, gy_key);
    gy_checked +:= 1;
end for;
// ---- PAIRED / CONIC PRESENTATIONS -------------------------------------------------------------
// Some bases are published, and stored by us, as a PAIR of equations rather than a single
// y^2=f(x) -- the curve is a complete intersection in a weighted projective space. Our model file
// stores these as <genus, "CRV", [Strings() | eq1, eq2]> with coordinates (x,y,z,s) and y of
// weight 2, i.e. P(1,2,1,1). Guo-Yang's own presentation may use different weights (21_2's second
// equation is degree 6, so its y has weight 3), which is fine -- IsIsomorphic compares the curves,
// not the embeddings.
Pours<x,y,z,s> := WeightedProjectiveSpace(Rationals(), [1,2,1,1]);

// <D, N, key, GY ambient weights, GY equations as a function of that space's coordinates>
gy_pairs := [*
    // X_0^21(2):  z^2 = -x^2-3  and  y^2 = -(3x-1)(3x+1)(x^2+7)(x^2+3)
    // homogenised with d; GY's y has weight 3 since its second equation has degree 6.  [genus 3]
    <21, 2, [Integers()|1], [1,3,1,1],
     func<a,b,c,d | [c^2 + a^2 + 3*d^2,
                     b^2 + (3*a-d)*(3*a+d)*(a^2+7*d^2)*(a^2+3*d^2)]>>,

    // X_0^57(1):  y^2 = (3s+1)(3s^3+11s^2+17s+1)  and  x^2 = -4s^2+2s-1               [genus 3]
    // GY's variables are (s,x,y), so here a=s, b=y, c=x, d=the homogenising variable; y has
    // weight 2 (its equation is degree 4). Involutions w_19(s,x,y)=(s,x,-y), w_57(s,x,y)=(s,-x,y).
    // ⚠ In the source this entry WRAPS across `\\`: the `x^2` line sits in a separate row of the
    // same \tabincell, so reading only the first math group silently yields a lone genus-1
    // y^2=quartic. Transcribed whole -- this is trap #1 of the four in the header.
    <57, 1, [Integers()|1], [1,2,1,1],
     func<a,b,c,d | [b^2 - (3*a+d)*(3*a^3+11*a^2*d+17*a*d^2+d^3),
                     c^2 + 4*a^2 - 2*a*d + d^2]>>
*];

for c in gy_pairs do
    gp_D, gp_N, gp_key, gp_w, gp_eqs := Explode(c);
    gp_mf := Sprintf("data/models/models_%o_%o.m", gp_D, gp_N);
    gp_models := eval (Read(gp_mf) cat "\nreturn models;");
    gp_ok, gp_entry := IsDefined(gp_models, gp_key);
    error if not gp_ok,
        Sprintf("X0^%o(%o): model file has no cover key %o", gp_D, gp_N, gp_key);
    gp_g, gp_tag, gp_strs := Explode(gp_entry[1]);
    error if Type(gp_tag) ne MonStgElt,
        Sprintf("X0^%o(%o) cover %o: expected a stored CRV pair, found a hyperelliptic model",
                gp_D, gp_N, gp_key);
    gp_Cours := Curve(Pours, [ eval ("return " cat str cat ";") : str in gp_strs ]);
    gp_Q<a,b,cc,d> := WeightedProjectiveSpace(Rationals(), gp_w);
    gp_Cgy := Curve(gp_Q, gp_eqs(a, b, cc, d));
    error if Genus(gp_Cours) ne Genus(gp_Cgy),
        Sprintf("X0^%o(%o) cover %o: our genus %o vs Guo-Yang's %o -- wrong object",
                gp_D, gp_N, gp_key, Genus(gp_Cours), Genus(gp_Cgy));
    error if not IsIsomorphic(gp_Cours, gp_Cgy),
        Sprintf("X0^%o(%o) cover %o: our CRV model is NOT isomorphic to Guo-Yang's published pair",
                gp_D, gp_N, gp_key);
    gy_checked +:= 1;
end for;

// ---------------------------------------------------------------------------------------------
// X_0^93(1): compared through the V_4 QUOTIENT DIAGRAM, not the full curve -- and the comparison
// determines a typo in the published table.
//
// GUO-YANG PRINT:   y^2 = (3s^3 - 7s^2 - 3t - 1)(3s^3 + s^2 - 3s - 9),   x^2 = -4s^2 - 6s - 9
// The `t` occurs nowhere else in that row, and the header above flags it as "almost certainly a
// typo". It is, and this test says WHICH typo: of the plausible repairs only `-3t -> -3s` yields a
// curve isomorphic to our genus-2 quotient. `-3s^2`, `-3`, and deleting the term all produce
// genus-2 curves that are NOT isomorphic to ours -- so the test discriminates rather than merely
// accommodating, which is what makes the conclusion worth anything. The three refuted repairs are
// checked below alongside the accepted one, precisely so that a future change cannot quietly turn
// this into a vacuous "some reading works" test. (The repo has produced three vacuous tests this
// way already; the guard is to count the comparisons actually made.)
//
// ⚠ WHY QUOTIENTS AND NOT THE FULL CURVE. The W={1} entry is a genus-5 CRV pair, and IsIsomorphic
// on those is the 10h+ regime measured at 26_3 -- not a CI cost. What IS checked: the two
// hyperelliptic quotients exactly, and the conic up to Q*/Q*^2. Three of the four cover keys are
// therefore pinned, and the fourth is their fibre product. Strong, but NOT a full-curve proof;
// models_93_1.m's header says the same. Closing it means an offline test like
// tests/_offline/GuoYangCurve_14_3.m.
gy93_P<gs> := PolynomialRing(Rationals());
gy93_B := 3*gs^3 + gs^2 - 3*gs - 9;                       // the second factor, unambiguous
gy93_C := -4*gs^2 - 6*gs - 9;                             // the conic
gy93_A := 3*gs^3 - 7*gs^2 - 3*gs - 1;                     // first factor, `-3t` read as `-3s`
gy93_models := eval (Read("data/models/models_93_1.m") cat "\nreturn models;");

// the accepted reading must match, on BOTH hyperelliptic quotients
for gy93_c in [* <[Integers()|1,93], gy93_A*gy93_B, "y-quotient">,
                 <[Integers()|1,3],  gy93_A*gy93_B*gy93_C, "product quotient"> *] do
    gy93_key, gy93_f, gy93_nm := Explode(gy93_c);
    gy93_ok, gy93_e := IsDefined(gy93_models, gy93_key);
    error if not gy93_ok, Sprintf("X0^93(1): model file has no cover key %o", gy93_key);
    gy93_ours := HyperellipticCurve(gy93_e[1][2]);
    gy93_gy   := HyperellipticCurve(gy93_f);
    error if Genus(gy93_ours) ne Genus(gy93_gy),
        Sprintf("X0^93(1) %o: our genus %o vs Guo-Yang's %o -- wrong object",
                gy93_nm, Genus(gy93_ours), Genus(gy93_gy));
    error if not IsIsomorphic(gy93_ours, gy93_gy),
        Sprintf("X0^93(1) %o (cover %o): NOT isomorphic to Guo-Yang's, reading `-3t` as `-3s`",
                gy93_nm, gy93_key);
    gy_checked +:= 1;
end for;

// the conic, up to squares: ours is -144s^2+36s-63, theirs -4s^2-6s-9; disc -34992 = -108*18^2.
gy93_oc := gy93_models[[Integers()|1,31]][1][2];
gy93_r  := Discriminant(gy93_oc) / Discriminant(gy93_C);
error if not IsSquare(Rationals()!gy93_r),
    Sprintf("X0^93(1) conic: disc ratio %o is not a square, so [1,31] is not Guo-Yang's conic",
            gy93_r);
gy_checked +:= 1;

// AND the three refuted repairs must STILL be refuted -- otherwise the typo is not determined.
gy93_ours93 := HyperellipticCurve(gy93_models[[Integers()|1,93]][1][2]);
gy93_alts := [* <"-3s^2", 3*gs^3 - 7*gs^2 - 3*gs^2 - 1>,
                <"dropped", 3*gs^3 - 7*gs^2 - 1>,
                <"-3",     3*gs^3 - 7*gs^2 - 3 - 1> *];
gy93_nref := 0;
for gy93_a in gy93_alts do
    gy93_nm, gy93_Aa := Explode(gy93_a);
    gy93_alt := HyperellipticCurve(gy93_Aa*gy93_B);
    error if Genus(gy93_alt) eq Genus(gy93_ours93) and IsIsomorphic(gy93_alt, gy93_ours93),
        Sprintf("X0^93(1): reading `-3t` as `%o` ALSO matches -- the typo is no longer determined, "
                * "and models_93_1.m's provenance claim must be weakened", gy93_nm);
    gy93_nref +:= 1;
end for;
error if gy93_nref ne 3,
    Sprintf("X0^93(1): expected 3 refuted readings, made %o comparisons", gy93_nref);

// gy_checked counts COMPARISONS, and 93_1 contributes three of them (two quotients + the conic)
// for one base -- so report both numbers rather than calling the total a base count.
printf " ok (%o comparison(s) over 10 base(s); 93_1 is quotient-level, %o alternative readings "
       * "refuted)\n", gy_checked, gy93_nref;
