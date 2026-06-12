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
    <34,  5,  {10,34},     "DN=170; #div(M=340)=12">,
    <34,  7,  {2,17},      "DN=238; #div(M=476)=12">,
    <74,  5,  {10,74},     "DN=370; #div(M=740)=12">,
    <10,  61, {10,122},    "DN=610; #div(M=1220)=12">,
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
    <119, 2,  {7,17},      "N=2: huge LP per Table 6">,
    <210, 19, {6,7,10,19}, "DN=3990: too large">
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
        crv_list, ws, keys := EquationsOfCovers(Xstar, curves);
        printf "  computed %o cover equations in %o s\n", #crv_list, Realtime()-t0;
        W := AllALsFromGens(gens, D*N);
        if not exists(k){k : k in keys | curves[k]`W eq W} then
            printf "  [D=%o N=%o W=<%o>] not among computed covers (keys); skipping\n", D, N, gens;
        else
            idx := Index(keys, k);
            C := crv_list[idx];
            check_group(C, gens, D, N);
        end if;
    catch e
        printf "  ERROR on (D,N)=(%o,%o): %o\n", D, N, e`Object;
    end try;
    printf "  ---- (D=%o,N=%o) done in %o s ----\n", D, N, Realtime()-t0;
end procedure;

curves := GetHyperellipticCandidates();
printf "Loaded %o candidate curves.\n", #curves;

if assigned idx then
    i := StringToInteger(idx);
    printf "Running single TABLE1 entry %o.\n", i;
    run_entry(TABLE1[i], curves);
else
    for entry in TABLE1 do
        run_entry(entry, curves);
    end for;
end if;

printf "\nDONE.\n";
exit;
