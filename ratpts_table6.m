// Check Table 6 (genus 0, unsure of X(Q)=empty) curves for rational points.
// For each (D,N) star curve, compute equations of immediate covers once,
// then for each target subgroup W (given by generator subscripts) check
// whether the genus-0 model has a rational point.

// Skip polymake LP instances that are too large to enumerate in reasonable time.
// 24*n is the bounding polytope coefficient; max solved n=499, D=51 N=2 hits n~78M.
LP_SIZE_CUTOFF := 10000;

AttachSpec("ShimuraQuotients.spec");
SetVerbose("ShimuraQuotients", 5);

curves := GetHyperellipticCandidates();
printf "Loaded %o candidate curves.\n", #curves;

// Each entry: <D, N, [ list of generator-sets, each a set of AL subscripts ] >
// CANDIDATES is set externally before loading, else defaults below.
if not assigned CANDIDATES then
    // squarefree-N only (method requires squarefree N); small-N batch.
    // DONE: D=10, N=7 — all three groups have rational points (genus 0 = conic ~ P^1):
    //   W=<{2,5}>:   27/16*x^2 + 47/64*x*z + y^2 + 5/64*z^2 = 0,  pt (-6/31 : 7/248 : 1)
    //   W=<{5,7}>:   27*x^2 - 22*x*z + y^2 - 5*z^2 = 0,            pt (-5/27 : 0 : 1)
    //   W=<{10,14}>: 27/64*x^2 + 5/64*x*z + y^2 = 0,               pt (-5/27 : 0 : 1)
    // NOTE: D=34, N=3 was attempted but failed — not enough CM points.
    // NOTE: D=26, N=5 was attempted but failed — "Could not find enough points".
    CANDIDATES := [*
        // --- priority cases to try next ---
        <21, 10, [ {2,15,21} ]>,
        // D*N = 102-255 (small-N batch)
        // N=2 cases deprioritized: D=51,N=2 was intractable (LP n~78M) and the
        // other N=2 cases likely share the same sparse-CM / huge-LP problem.
        <55, 2, [ {2,55} ]>,
        <87, 2, [ {2,87} ]>,
        <95, 2, [ {10,38} ]>,
        <111, 2, [ {2,111} ]>,
        <51, 5, [ {5,51} ]>,
        <35, 6, [ {2,15,21} ]>,
        // D*N = 210
        <15, 14, [ {6,7,10} ]>,
        <14, 15, [ {3,10,14}, {5,6,14} ]>,
        <10, 21, [ {3,5,7}, {3,10,14} ]>,
        <6,  35, [ {5,6,14}, {6,7,10} ]>,
        // D*N = 330
        <33, 10, [ {2,15,33} ]>,
        <22, 15, [ {6,10,11} ]>,
        <15, 22, [ {3,10,22} ]>,
        <10, 33, [ {2,15,33} ]>,
        <6,  55, [ {3,10,22}, {6,10,11} ]>,
        // D*N = 462
        <21, 22, [ {6,7,22} ]>,
        <22, 21, [ {3,14,22} ]>,
        <6,  77, [ {2,21,33}, {6,11,14} ]>,
        // D*N = 510
        <15, 34, [ {2,15,17} ]>,
        // D*N = 546
        <14, 39, [ {3,13,14} ]>,
        <26, 21, [ {3,7,26} ]>,
        // D*N = 690
        <10, 69, [ {3,5,46} ]>,
        // D*N = 770
        <10, 77, [ {7,10,11} ]>,
        // D*N = 2730 (likely too large; LP cutoff will bail quickly)
        <390, 7, [ {3,7,10,13} ]>
    *];
end if;

// Ordered, indexable view of the Table 6 groups (D*N ascending), for run_table.sh.
TABLE6 := [ x : x in CANDIDATES ];  // convert List to SeqEnum for Sort
Sort(~TABLE6, func< a, b | (a[1]*a[2]) ne (b[1]*b[2]) select (a[1]*a[2])-(b[1]*b[2])
                           else (a[2] ne b[2] select a[2]-b[2] else a[1]-b[1]) >);

if assigned idx then
    i := StringToInteger(idx);
    printf "Running single TABLE6 group #%o of %o.\n", i, #TABLE6;
    CANDIDATES := [* TABLE6[i] *];   // restrict the main loop to this one group
end if;

procedure check_group(C, gens, D, N)
    g := Genus(C);
    desc := Sprintf("D=%o N=%o W=<%o>", D, N, gens);
    if g ne 0 then
        printf "  [%o] WARNING genus = %o (expected 0)\n", desc, g;
        return;
    end if;
    f := HyperellipticPolynomials(C);
    if Degree(f) le 1 then
        printf "  [%o] model is P^1 => HAS rational point\n", desc;
        return;
    end if;
    con := Conic(C);
    has, pt := HasRationalPoint(con);
    if has then
        printf "  [%o] HAS rational point: %o   (conic: %o)\n", desc, pt, con;
    else
        printf "  [%o] NO rational point (conic pointless): %o\n", desc, con;
    end if;
end procedure;

// polymake LP dimension guard. The Borcherds-form step enumerates lattice
// points of a polytope of dimension #Divisors(M), M = 4*D0 = 4*(D*N)/2^v2(D)
// (the level of the {oo}-weakly-holomorphic forms ring). Cost is driven by
// #Divisors(M), NOT by LP_SIZE_CUTOFF (which bounds the pole order n only):
// the OOM-killed M=420 case had n=145 << 10000. Empirically (completed vs
// Killed:9 polymake artifacts): #div<=12 completes reliably (proven to M=1212,
// n=499); #div=16-20 only at small n; #div>=24 (3 odd prime factors of M) OOMs
// even at its minimum forced n. So skip #div(M) >= DIV_CUTOFF up front -- this
// turns the kills into clean skips; it does NOT make new cases tractable.
DIV_CUTOFF := 24;
polymake_level := func< D, N | 4 * ((D*N) div 2^Valuation(D, 2)) >;  // = M

for entry in CANDIDATES do
    D := entry[1]; N := entry[2]; gensets := entry[3];
    printf "\n==== D=%o N=%o (D*N=%o) ====\n", D, N, D*N;
    if not IsSquarefree(N) then
        printf "  N=%o is not squarefree; method N/A; skipping\n", N;
        continue;
    end if;
    M := polymake_level(D, N);
    ndiv := #Divisors(M);
    if ndiv ge DIV_CUTOFF then
        printf "  polymake level M=%o has #div=%o >= %o; OOM-doomed, skipping\n",
            M, ndiv, DIV_CUTOFF;
        continue;
    end if;
    t0 := Realtime();
    if not exists(Xstar){X : X in curves | X`D eq D and X`N eq N and IsStarCurve(X)} then
        printf "  no star curve found for (D,N)=(%o,%o); skipping\n", D, N;
        continue;
    end if;
    try
        crv_list, ws, keys := EquationsOfCovers(Xstar, curves);
        printf "  computed %o cover equations in %o s\n", #crv_list, Realtime()-t0;
        for gens in gensets do
            W := AllALsFromGens(gens, D*N);
            if not exists(k){k : k in keys | curves[k]`W eq W} then
                printf "  [D=%o N=%o W=<%o>] not among computed covers (keys); skipping\n", D, N, gens;
                continue;
            end if;
            idx := Index(keys, k);
            C := crv_list[idx];
            check_group(C, gens, D, N);
        end for;
    catch e
        printf "  ERROR on (D,N)=(%o,%o): %o\n", D, N, e`Object;
    end try;
    printf "  ---- group (D=%o,N=%o) done in %o s ----\n", D, N, Realtime()-t0;
end for;

printf "\nDONE.\n";
exit;
