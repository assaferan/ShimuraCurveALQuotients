// Check Table 6 (genus 0, unsure of X(Q)=empty) curves for rational points.
// For each (D,N) star curve, compute equations of immediate covers once,
// then for each target subgroup W (given by generator subscripts) check
// whether the genus-0 model has a rational point.

AttachSpec("ShimuraQuotients.spec");
SetVerbose("ShimuraQuotients", 1);

curves := GetHyperellipticCandidates();
printf "Loaded %o candidate curves.\n", #curves;

// Each entry: <D, N, [ list of generator-sets, each a set of AL subscripts ] >
// CANDIDATES is set externally before loading, else defaults below.
if not assigned CANDIDATES then
    // squarefree-N only (method requires squarefree N); small-N batch.
    CANDIDATES := [*
        <51, 2, [ {6,34} ]>,
        <55, 2, [ {2,55} ]>,
        <87, 2, [ {2,87} ]>,
        <95, 2, [ {10,38} ]>,
        <111, 2, [ {2,111} ]>,
        <26, 5, [ {10,26} ]>,
        <51, 5, [ {5,51} ]>,
        <35, 6, [ {2,15,21} ]>
    *];
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

for entry in CANDIDATES do
    D := entry[1]; N := entry[2]; gensets := entry[3];
    printf "\n==== D=%o N=%o (D*N=%o) ====\n", D, N, D*N;
    if not IsSquarefree(N) then
        printf "  N=%o is not squarefree; method N/A; skipping\n", N;
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
