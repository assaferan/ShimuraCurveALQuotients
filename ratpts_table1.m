// Check Table 1 (genus 1, authors unsure whether X(Q)=empty) curves for
// rational points. For each (D,N) star curve, compute equations of immediate
// covers once, then for the target subgroup W (given by generator subscripts)
// pull out the genus-1 model C and test IsEllipticCurve(C). If C is an
// elliptic curve it HAS a rational point (resolving X(Q) != empty); the
// returned E / the hyperelliptic involution = image of Fricke w_{DN}.
//
// This is the genus-1 analog of ratpts_table6.m (which handled genus-0 conics).
//
// Usage:
//   magma ratpts_table1.m                 // run the full ordered TABLE1
//   magma idx:=2 ratpts_table1.m          // run only TABLE1[2]
//
// Resource notes / running log live in the RESULTS LOG block at the bottom and
// are appended to ratpts_table1_output.txt as runs complete. Once a (D,N,W) has
// a recorded MODEL or a recorded reason-it-fails, do NOT re-run it.

//Errors:
//  ERROR on (D,N)=(34,5): Could not find enough points, sorry!
AttachSpec("ShimuraQuotients.spec");
SetVerbose("ShimuraQuotients", 1);

// ---------------------------------------------------------------------------
// Table 1 entries as <D, N, gens, note>, ordered by expected feasibility
// (roughly D*N ascending, N>=5 / moderate-D first; tiny-N=1,2 and very large
// D*N deprioritized to the end because they hit the sparse-CM / huge-LP wall
// seen for the N=2,3 cases in Table 6).
//   gens = set of AL subscripts generating W  (e.g. <w10,w42>  ->  {10,42})
// ---------------------------------------------------------------------------
TABLE1 := [*
    // ---- TRACTABLE batch: #div(M)<=12 (M=4*p*q), the only rows that complete ----
    <34,  5,  {10,34},     "DN=170; #div(M=340)=12 | 2026-06-15 SKIP insufficient CM (need 7, have 4)">,
    <34,  7,  {2,17},      "DN=238; #div(M=476)=12 | 2026-06-15 SKIP insufficient CM (need 7, have 5)">,
    <74,  5,  {10,74},     "DN=370; #div(M=740)=12 | 2026-06-15 SKIP insufficient CM (need 7, have 4)">,
    <10,  61, {10,122},    "DN=610; #div(M=1220)=12 | 2026-06-15 TIMEOUT@900s (reached polymake (1220,513,0); UNRESOLVED, retry w/ bigger budget)">,
    // ---- #div(M)>=24: OOM wall, skipped by DIV_CUTOFF (kept for record) ----
    <6,   35, {10,42},     "DN=210; #div(M=420)=24 OOM">,
    <10,  21, {5,21},      "DN=210; #div24 OOM">,
    <15,  14, {2,5,7},     "DN=210; #div(M=840)=32 OOM">,
    <10,  33, {2,55},      "DN=330; #div24 OOM">,
    <10,  39, {3,5,13},    "DN=390; #div24 OOM">,
    <10,  51, {2,15,51},   "DN=510; #div24 OOM">,
    <6,   85, {2,15,51},   "DN=510; #div24 OOM">,
    <21,  26, {2,21,39},   "DN=546; #div(M=2184)=32 OOM">,
    <6,   115,{5,6,46},    "DN=690; #div24 OOM">,
    <14,  57, {2,21,57},   "DN=798; #div24 OOM">,
    <10,  93, {3,10,62},   "DN=930; #div24 OOM">,
    <6,   161,{2,21,69},   "DN=966; #div24 OOM">,
    // ---- deprioritized: N=1/2 (sparse CM, huge LP) and/or very large D*N ----
    <210, 1,  {7,15},      "DN=210 but N=1: sparse CM, likely intractable">,
    <330, 1,  {2,33},      "N=1">,
    <330, 1,  {3,10},      "N=1, second W">,
    <462, 1,  {11,14},     "N=1">,
    <798, 1,  {2,3,19},    "N=1">,
    <1230,1,  {3,10,82},   "N=1, large D">,
    <1722,1,  {6,14,41},   "N=1, large D">,
    <119, 2,  {7,17},      "N=2: huge LP per Table 6 | 2026-06-15 SKIP insufficient CM (need 7, have 1)">,
    <210, 19, {6,7,10,19}, "DN=3990: too large">
*];

// ---------------------------------------------------------------------------
// Table 7 (from OanaFreddy.pdf): 54 genus-1 curves X = X_0^D(N)/W with Jac of
// positive rank, authors unsure if X(Q)=empty. SAME check as Table 1 -- only
// (D,N,gens) is needed (the model is computed, not read), so no f column.
// Ordered tractable-first (#div(M)<=12, squarefree N), then the OOM-wall rest.
// Most rows are #div(M)>=24 and will be skipped by the DIV_CUTOFF guard.
// ---------------------------------------------------------------------------
TABLE7 := [*
    <134, 3, {3,134}, "T7 DN=402; M=804 #div12 | 2026-06-15 SKIP insufficient CM (need 7, have 1)">,
    <122, 7, {2,7,61}, "T7 DN=854; M=1708 #div12 | 2026-06-15 SKIP: W not an immediate cover of X* (not in candidate data)">,
    <34, 29, {2,17,29}, "T7 DN=986; M=1972 #div12 | 2026-06-15 SKIP: W not an immediate cover of X* (not in candidate data)">,
    <26, 9, {9,26}, "T7 DN=234; M=468 #div18 NOT-sqfree">,
    <95, 3, {3,95}, "T7 DN=285; M=1140 #div24">,
    <159, 2, {2,159}, "T7 DN=318; M=1272 #div16 | 2026-06-15 SKIP insufficient CM (need 7, have 0)">,
    <14, 25, {14,25}, "T7 DN=350; M=700 #div18 NOT-sqfree">,
    <39, 10, {2,5,39}, "T7 DN=390; M=1560 #div32">,
    <15, 26, {2,13,15}, "T7 DN=390; M=1560 #div32">,
    <215, 2, {2,5,43}, "T7 DN=430; M=1720 #div16 | 2026-06-15 SKIP: W not an immediate cover of X* (not in candidate data)">,
    <77, 6, {2,3,7,11}, "T7 DN=462; M=1848 #div32">,
    <33, 14, {2,7,33}, "T7 DN=462; M=1848 #div32">,
    <10, 49, {5,98}, "T7 DN=490; M=980 #div18 NOT-sqfree">,
    <85, 6, {2,3,5,17}, "T7 DN=510; M=2040 #div32">,
    <51, 10, {2,5,51}, "T7 DN=510; M=2040 #div32">,
    <95, 6, {2,3,5,19}, "T7 DN=570; M=2280 #div32">,
    <15, 38, {2,3,5,19}, "T7 DN=570; M=2280 #div32">,
    <22, 35, {5,7,22}, "T7 DN=770; M=1540 #div24">,
    <35, 26, {2,5,7,13}, "T7 DN=910; M=3640 #div32">,
    <6, 155, {2,5,93}, "T7 DN=930; M=1860 #div24">,
    <6, 155, {5,6,31}, "T7 DN=930; M=1860 #div24">,
    <6, 169, {2,3,169}, "T7 DN=1014; M=2028 #div18 NOT-sqfree">,
    <22, 51, {2,3,11,17}, "T7 DN=1122; M=2244 #div24">,
    <21, 55, {3,5,7,11}, "T7 DN=1155; M=4620 #div48">,
    <6, 203, {3,7,58}, "T7 DN=1218; M=2436 #div24">,
    <15, 82, {2,3,5,41}, "T7 DN=1230; M=4920 #div32">,
    <14, 95, {2,5,7,19}, "T7 DN=1330; M=2660 #div24">,
    <10, 141, {2,3,5,47}, "T7 DN=1410; M=2820 #div24">,
    <38, 39, {2,3,13,19}, "T7 DN=1482; M=2964 #div24">,
    <10, 159, {2,3,5,53}, "T7 DN=1590; M=3180 #div24">,
    <10, 161, {2,5,7,23}, "T7 DN=1610; M=3220 #div24">,
    <6, 287, {2,3,7,41}, "T7 DN=1722; M=3444 #div24">,
    <1155, 2, {2,3,5,7,11}, "T7 DN=2310; M=9240 #div64">,
    <770, 3, {2,3,5,77}, "T7 DN=2310; M=4620 #div48">,
    <770, 3, {2,3,7,55}, "T7 DN=2310; M=4620 #div48">,
    <770, 3, {3,5,11,14}, "T7 DN=2310; M=4620 #div48">,
    <770, 3, {3,7,10,11}, "T7 DN=2310; M=4620 #div48">,
    <770, 3, {3,10,14,22}, "T7 DN=2310; M=4620 #div48">,
    <10, 231, {2,3,5,7,11}, "T7 DN=2310; M=4620 #div48">,
    <546, 5, {2,3,5,91}, "T7 DN=2730; M=5460 #div48">,
    <546, 5, {2,5,13,21}, "T7 DN=2730; M=5460 #div48">,
    <546, 5, {3,5,7,26}, "T7 DN=2730; M=5460 #div48">,
    <546, 5, {3,5,13,14}, "T7 DN=2730; M=5460 #div48">,
    <546, 5, {5,6,14,26}, "T7 DN=2730; M=5460 #div48">,
    <390, 7, {2,3,7,65}, "T7 DN=2730; M=5460 #div48">,
    <390, 7, {2,5,7,39}, "T7 DN=2730; M=5460 #div48">,
    <390, 7, {2,7,13,15}, "T7 DN=2730; M=5460 #div48">,
    <6, 455, {2,3,5,7,13}, "T7 DN=2730; M=5460 #div48">,
    <1190, 3, {2,3,5,7,17}, "T7 DN=3570; M=7140 #div48">,
    <570, 7, {2,3,7,95}, "T7 DN=3990; M=7980 #div48">,
    <1430, 3, {2,3,5,11,13}, "T7 DN=4290; M=8580 #div48">,
    <858, 5, {2,3,5,11,13}, "T7 DN=4290; M=8580 #div48">,
    <510, 11, {2,3,5,11,17}, "T7 DN=5610; M=11220 #div48">,
    <546, 11, {2,3,7,11,13}, "T7 DN=6006; M=12012 #div48">
*];

procedure check_group(C, gens, D, N)
    desc := Sprintf("D=%o N=%o W=<%o>", D, N, gens);
    g := Genus(C);
    if g ne 1 then
        printf "  [%o] WARNING genus = %o (expected 1)\n", desc, g;
        return;
    end if;
    f := HyperellipticPolynomials(C);
    printf "  [%o] model y^2 = %o\n", desc, f;
    is_ell, E := IsEllipticCurve(C);
    if is_ell then
        Em := MinimalModel(E);
        printf "  [%o] IS elliptic => HAS rational point.\n", desc;
        printf "  [%o]   E (min model): %o  conductor=%o  rank(an?)=...\n", desc, Em, Conductor(Em);
        printf "  [%o]   aInvariants=%o\n", desc, aInvariants(Em);
    else
        printf "  [%o] genus 1, IsEllipticCurve found NO rational point (model only).\n", desc;
    end if;
end procedure;

// polymake LP dimension guard (see ratpts_table6.m / HANDOFF_table6.md for the
// full #div(M) tractability law). The Borcherds-form step enumerates lattice
// points of a polytope of dimension #Divisors(M), M = 4*(D*N)/2^v2(D). Cost is
// driven by #div(M), NOT by the pole order n (so LP_SIZE_CUTOFF misses it: the
// killed M=420 case had n=145). #div<=12 completes; #div>=24 (3 odd primes in M)
// OOMs even at minimum forced n. Skip those cleanly instead of `Killed: 9`.
DIV_CUTOFF := 24;
polymake_level := func< D, N | 4 * ((D*N) div 2^Valuation(D, 2)) >;  // = M

procedure run_entry(entry, curves)
    D := entry[1]; N := entry[2]; gens := entry[3]; note := entry[4];
    printf "\n==== D=%o N=%o (D*N=%o) W=<%o>  [%o] ====\n", D, N, D*N, gens, note;
    if not IsSquarefree(N) then
        printf "  N=%o is not squarefree; method N/A; skipping\n", N;
        return;
    end if;
    M := polymake_level(D, N);
    ndiv := #Divisors(M);
    if ndiv ge DIV_CUTOFF then
        printf "  polymake level M=%o has #div=%o >= %o; OOM-doomed, skipping\n",
            M, ndiv, DIV_CUTOFF;
        return;
    end if;
    t0 := Realtime();
    if not exists(Xstar){X : X in curves | X`D eq D and X`N eq N and IsStarCurve(X)} then
        printf "  no star curve found for (D,N)=(%o,%o); skipping\n", D, N;
        return;
    end if;
    try
        W := AllALsFromGens(gens, D*N);
        tgts := { W };
        // If the target W isn't an immediate cover of X* in the candidate data, skip
        // cleanly (don't let the Targets require throw) -- nothing to compute.
        if not exists(k0){i : i in Xstar`CoveredBy | curves[i]`W eq W} then
            printf "  [D=%o N=%o W=<%o>] not an immediate cover of X* (not in candidate data); skipping\n", D, N, gens;
        else
            // Cheap predict-and-skip on CM-point count, restricted to this target.
            enough, need, have := EnoughCMPointsForTargets(Xstar, curves, tgts);
            if not enough then
                printf "  insufficient CM points (need=%o, have=%o); skipping before Borcherds work\n", need, have;
            else
                // Targets-restricted: build only this cover (lower num_vals, rescues a
                // genus-1 target inflated by genus-2 siblings; speedup otherwise).
                crv_list, ws, keys := EquationsOfCovers(Xstar, curves : Targets := tgts);
                printf "  computed %o cover equations (targets-restricted) in %o s\n", #crv_list, Realtime()-t0;
                if not exists(k){k : k in keys | curves[k]`W eq W} then
                    printf "  [D=%o N=%o W=<%o>] not among computed covers (keys); skipping\n", D, N, gens;
                else
                    idx := Index(keys, k);
                    C := crv_list[idx];
                    check_group(C, gens, D, N);
                end if;
            end if;
        end if;
    catch e
        printf "  ERROR on (D,N)=(%o,%o): %o\n", D, N, e`Object;
    end try;
    printf "  ---- (D=%o,N=%o) done in %o s ----\n", D, N, Realtime()-t0;
end procedure;

SetVerbose("ShimuraQuotients",5);
curves := GetHyperellipticCandidates();
printf "Loaded %o candidate curves.\n", #curves;

// Table 1 and Table 7 use the identical check; run them as one ordered list.
ALL := TABLE1 cat TABLE7;
printf "Loaded %o Table 1 + %o Table 7 = %o entries.\n", #TABLE1, #TABLE7, #ALL;

if assigned idx then
    i := StringToInteger(idx);
    printf "Running single entry %o (of %o).\n", i, #ALL;
    run_entry(ALL[i], curves);
elif assigned table and StringToInteger(table) eq 7 then
    printf "Running TABLE7 only (%o entries).\n", #TABLE7;
    for entry in TABLE7 do
        run_entry(entry, curves);
    end for;
else
    for entry in ALL do
        run_entry(entry, curves);
    end for;
end if;

printf "\nDONE.\n";
exit;
