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
// z, and 93_1 mixes two variables -- almost certainly a typo in the paper. Do not re-automate this.
//
// ⚠ AND THE PUBLISHED EQUATION IS NOT ALWAYS OUR W={1} KEY. `X0_26_1.m` compares a degree-3
// s-equation against Guo-Yang's degree-6 x-equation; they differ by s = x^2. So each case below
// records WHICH cover key it is comparing, established per base rather than assumed.
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
    <51, 1, [Integers()|1], -(x^2+3)*(243*x^6+235*x^4-31*x^2+1)>
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
printf " ok (%o base(s))\n", gy_checked;
