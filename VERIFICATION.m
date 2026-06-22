// ============================================================================
// VERIFICATION.m
//
// Reconstructs, from scratch, the genus-2 quotient models published in the two
// paper tables ("modelsgen2" and "modelsbielliptic") and checks that each is
// ISOMORPHIC (over Q) to the printed model  y^2 + h(x) y = f(x).
//
// Each row is a degree-2 cover of a star curve X(D,N)*. The published equations
// use different hauptmoduls per row (see the paper's remark), so the literal
// polynomials are NOT expected to match -- only the isomorphism class is.
//
// Two cover depths occur and are handled like eqns_genus2_bielliptic.m:
//   * W an IMMEDIATE cover of X*  -> EquationsOfCovers(Xstar, curves : Targets)
//   * W a DEEPER cover            -> AllEquationsAboveCovers (build upward)
// (The genus-2 non-bielliptic rows are all immediate; the two (6,17) bielliptic
//  rows sit two rungs above X* and need the climb.)
//
// USAGE:
//   magma VERIFICATION.m            // verify all rows
//   magma row:=7 VERIFICATION.m     // verify only row #7
//   magma lo:=1 hi:=8 VERIFICATION.m// verify rows #1..#8
//
// Uses the cached polymake solutions in polymake/, so polymake/libnormaliz is
// not required to run this.
//
// RUNTIME: ~3-9 min per row (dominated by the per-star-curve Borcherds-form
// build), so the full 17-row sweep takes on the order of 1-2 hours. Use the
// row:= / lo:= / hi:= selectors to verify a subset.
// ============================================================================

AttachSpec("ShimuraQuotients.spec");
SetVerbose("ShimuraQuotients", 1);

R<x> := PolynomialRing(Rationals());

// Each row: < D, N, gens, h, f, note >  modelling  y^2 + h*y = f.
// `gens` is the set of Atkin-Lehner subscripts generating W.
// NOTE: row (6,37) is published as <w_6, w_34>, but w_34 is not a valid AL
// involution of X^6(37) (34 is not a Hall divisor of 222); read as a typo for
// w_37 here -- see `note`.
ROWS := [*
    < 6,  29, {3,29},   -x^2-x,            -144*x^5-117*x^4+41*x^3+21*x^2-6*x,                                          "" >,
    < 14, 13, {2,13},   2*x^3-3*x^2+x,     -44*x^6-626*x^5-3296*x^4-8298*x^3-10950*x^2-7306*x-1950,                    "" >,
    < 14, 13, {7,13},   2*x^3-3*x^2+x,     -1716*x^6+1816*x^5-117*x^4-291*x^3+11*x^2+18*x+2,                           "" >,
    < 10, 17, {10,34},  x^2-1,             -4*x^5+11*x^4-33*x^3+21*x^2-15*x-44,                                         "" >,
    < 14, 29, {7,29},   -x^2+x-2,          1792*x^5-6487*x^4+9347*x^3-6701*x^2+2390*x-340,                             "" >,
    < 14, 29, {7,58},   -x-1,              -224*x^5+365*x^4-194*x^3+47*x^2-6*x,                                         "" >,
    < 6,  37, {6,37},   -x^2-x,            972*x^5-2411*x^4+2244*x^3-929*x^2+144*x,                                     "table prints w_34 (invalid); read as w_37" >,
    < 6,  41, {2,3},    -x^2-x,            -2304*x^6-3862*x^5-2270*x^4-687*x^3-230*x^2-76*x-10,                        "" >,
    < 6,  41, {3,41},   -x^2-x,            -288*x^5-303*x^4-95*x^3-27*x^2-12*x-2,                                       "" >,
    < 6,  43, {2,3},    -x^2-1,            -15552*x^6+24614*x^5-16682*x^4+5999*x^3-1182*x^2+119*x-5,                   "" >,
    < 10, 53, {5,106},  -x^2-x,            20*x^5-25*x^4+35*x^3-16*x^2+9*x+1,                                           "" >,
    < 6,  59, {2,59},   -x^2-x,            6*x^5+16*x^4-23*x^3-26*x^2+84*x-40,                                          "" >,
    < 6,  61, {2,3},    x^2-1,             -559872*x^6+2603694*x^5-5041398*x^4+5202037*x^3-3016984*x^2+932433*x-119975, "" >,
    < 106, 1, {2},      -2*x^3+3*x^2-3*x,  -20*x^6+164*x^5-540*x^4+872*x^3-700*x^2+250*x-32,                           "" >,
    < 118, 1, {2},      -2*x^3-x^2,        -688676*x^6+1519312*x^5-1397890*x^4+686641*x^3-189920*x^2+28048*x-1728,     "" >,
    // ---- modelsbielliptic table (genus 2, bielliptic AL involution) ----
    < 6,  17, {3},      -x^3-x,            -13*x^6-61*x^4+33*x^2+36,                                                    "bielliptic" >,
    < 6,  17, {6},      0,                 -9*x^6+27*x^5-94*x^4+143*x^3-446*x^2+379*x-321,                              "bielliptic" >
*];

// Select which rows to run.
if assigned row then
    SEL := [StringToInteger(row)];
elif assigned lo then
    hi_ := assigned hi select StringToInteger(hi) else #ROWS;
    SEL := [StringToInteger(lo)..hi_];
else
    SEL := [1..#ROWS];
end if;

curves := GetHyperellipticCandidates();
printf "Loaded %o candidate curves. Verifying %o row(s).\n\n", #curves, #SEL;

// Compute the model curve C for a single row, or return false + reason.
// Returns < ok, C, reason >.
function compute_model(D, N, gens, curves)
    if not exists(Xstar){X : X in curves | X`D eq D and X`N eq N and IsStarCurve(X)} then
        return false, _, "no star curve found";
    end if;
    W := AllALsFromGens(gens, D*N);
    imm_W := { curves[i]`W : i in Xstar`CoveredBy };
    if W in imm_W then
        // Immediate cover: target-restricted EquationsOfCovers.
        crv_list, ws, keys := EquationsOfCovers(Xstar, curves : Targets := {W});
        if not exists(k){k : k in keys | curves[k]`W eq W} then
            return false, _, "immediate cover not among computed keys";
        end if;
        return true, crv_list[Index(keys, k)], "immediate";
    else
        // Deeper cover: build upward.
        all_eqns, all_ws := AllEquationsAboveCovers(D, N, curves);
        if not exists(k){k : k in Keys(all_eqns) |
                         curves[k]`D eq D and curves[k]`N eq N and curves[k]`W eq W} then
            return false, _, "no curve id with this W in all_eqns";
        end if;
        if #Keys(all_eqns[k]) eq 0 then
            return false, _, "cover not reachable by climbing (no base eqn)";
        end if;
        return true, all_eqns[k][Representative(Keys(all_eqns[k]))], "climb";
    end if;
end function;

n_match := 0;
verdicts := [];   // < idx, desc, tag >
models := [* *];  // < idx, desc, C > for every row whose model was computed

for i in SEL do
    rrow := ROWS[i];
    D := rrow[1]; N := rrow[2]; gens := rrow[3]; h := rrow[4]; f := rrow[5]; note := rrow[6];
    desc := Sprintf("#%o X(%o,%o)/<%o>", i, D, N, gens);
    if note ne "" then desc cat:= Sprintf(" [%o]", note); end if;
    Ctab := HyperellipticCurve(f, h);   // y^2 + h*y = f
    printf "==== %o ====\n", desc;
    t0 := Realtime();
    try
        ok, C, reason := compute_model(D, N, gens, curves);
        if ok then Append(~models, <i, desc, C>); end if;
        if not ok then
            printf "  COULD-NOT-COMPUTE: %o\n", reason;
            Append(~verdicts, <i, desc, "NO-MODEL">);
        elif IsIsomorphic(C, Ctab) then
            printf "  MATCH (%o) in %o s\n", reason, Realtime()-t0;
            n_match +:= 1;
            Append(~verdicts, <i, desc, "MATCH">);
        else
            same := G2Invariants(C) eq G2Invariants(Ctab);
            printf "  MISMATCH (%o); G2Invariants %o\n", reason,
                   same select "agree (Qbar-isomorphic)" else "differ";
            printf "    computed: y^2 + (%o)*y = %o\n", HyperellipticPolynomials(C);
            Append(~verdicts, <i, desc, same select "MISMATCH(geom-iso)" else "MISMATCH">);
        end if;
    catch e
        printf "  ERROR: %o\n", e`Object;
        Append(~verdicts, <i, desc, "ERROR">);
    end try;
    printf "\n";
end for;

printf "================ SUMMARY ================\n";
for v in verdicts do
    printf "  %-7o %o\n", v[3], v[2];
end for;
printf "----------------------------------------\n";
printf "  %o / %o matched\n", n_match, #SEL;
printf "========================================\n";

// ---- Geometric automorphism groups (brief check over the computed models) ----
// |Aut| is 2 generically (just the hyperelliptic involution) and 4 for the
// bielliptic curves: X(14,29)/<w7,w29> and both D=6,N=17 rows.
printf "\n========= GEOMETRIC AUTOMORPHISM GROUPS =========\n";
aut_all_ok := true;
for m in models do
    rrow := ROWS[m[1]];
    D := rrow[1]; N := rrow[2]; gens := rrow[3];
    exp_aut := ((D eq 14 and N eq 29 and gens eq {7,29}) or (D eq 6 and N eq 17))
               select 4 else 2;
    naut := #GeometricAutomorphismGroup(m[3]);
    ok := naut eq exp_aut;
    aut_all_ok := aut_all_ok and ok;
    printf "  %-6o %o : |Aut| = %o (expected %o)\n", ok select "OK" else "WRONG",
           m[2], naut, exp_aut;
end for;
printf "-------------------------------------------------\n";
assert aut_all_ok;
printf "  all geometric automorphism group sizes as expected.\n";
printf "=================================================\n";

exit;
