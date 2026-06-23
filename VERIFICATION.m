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
// Rows are GROUPED per star curve (D,N): all of a curve's immediate targets are
// solved in ONE EquationsOfCovers call, and its deep targets in ONE
// AllEquationsAboveCovers call. This runs the expensive Borcherds-form build once
// per star curve instead of once per row.
//
// Each row's isomorphism to the table model is checked inline, right after its
// star curve is built (verdicts appear incrementally). Then a SEPARATE final pass
// checks ALL computed models for biellipticity via the geometric automorphism
// group: |Aut| = 2 (non-bielliptic) generically, and |Aut| = 4 (bielliptic) for
// X(14,29)/<w7,w29> and both D=6,N=17 rows.
//
// USAGE:
//   magma VERIFICATION.m            // verify all rows
//   magma row:=7 VERIFICATION.m     // verify only row #7
//   magma lo:=1 hi:=8 VERIFICATION.m// verify rows #1..#8
//
// Uses the cached polymake solutions in polymake/, so polymake/libnormaliz is
// not required to run this.
//
// RUNTIME: the cost is one Borcherds-form build per distinct star curve (the
// 17 rows cover 13 star curves), ranging from ~2 min to ~25 min each (heaviest:
// levels M=812 and M=1060), so the full grouped sweep takes on the order of
// 2-3 hours. Use the row:= / lo:= / hi:= selectors to verify a subset.
// ============================================================================

AttachSpec("ShimuraQuotients.spec");
SetVerbose("ShimuraQuotients", 1);

R<x> := PolynomialRing(Rationals());

// Each row: < D, N, gens, h, f, note >  modelling  y^2 + h*y = f.
// `gens` is the set of Atkin-Lehner subscripts generating W.
// NOTE: row (6,37) is published as <w_6, w_34>, but w_34 is not a valid AL
// involution of X^6(37) (34 is not a Hall divisor of 222); read as a typo for
// w_74 here -- <w_6,w_74> is the genus-2 cover, whereas <w_6,w_37> is genus 0.
ROWS := [*
    < 6,  29, {3,29},   -x^2-x,            -144*x^5-117*x^4+41*x^3+21*x^2-6*x,                                          "" >,
    < 14, 13, {2,13},   2*x^3-3*x^2+x,     -44*x^6-626*x^5-3296*x^4-8298*x^3-10950*x^2-7306*x-1950,                    "" >,
    < 14, 13, {7,13},   2*x^3-3*x^2+x,     -1716*x^6+1816*x^5-117*x^4-291*x^3+11*x^2+18*x+2,                           "" >,
    < 10, 17, {10,34},  x^2-1,             -4*x^5+11*x^4-33*x^3+21*x^2-15*x-44,                                         "" >,
    < 14, 29, {7,29},   -x^2+x-2,          1792*x^5-6487*x^4+9347*x^3-6701*x^2+2390*x-340,                             "" >,
    < 14, 29, {7,58},   -x-1,              -224*x^5+365*x^4-194*x^3+47*x^2-6*x,                                         "" >,
    < 6,  37, {6,74},   -x^2-x,            972*x^5-2411*x^4+2244*x^3-929*x^2+144*x,                                     "table prints w_34 (invalid); read as w_74 (genus 2; w_37 gives genus 0)" >,
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
// Optional: skip:=11 (or skip:="11,7") removes those rows from the selection.
// Row 11, X(10,53) (M=1060), is impractically slow in its [4/6] step (>13h),
// so skip:=11 lets the other rows complete.
if assigned skip then
    skipset := { StringToInteger(s) : s in Split(skip, ",") };
    SEL := [i : i in SEL | i notin skipset];
end if;

// Short label for a row (used in all reporting).
function row_desc(i)
    rr := ROWS[i];
    d := Sprintf("#%o X(%o,%o)/<%o>", i, rr[1], rr[2], rr[3]);
    if rr[6] ne "" then d cat:= Sprintf(" [%o]", rr[6]); end if;
    return d;
end function;

// Expected geometric automorphism group size: 2 generically (just the
// hyperelliptic involution), 4 for the bielliptic curves -- X(14,29)/<w7,w29>
// and both D=6,N=17 rows (which carry an extra, bielliptic involution).
function expected_aut(i)
    rr := ROWS[i]; D := rr[1]; N := rr[2]; gens := rr[3];
    return ((D eq 14 and N eq 29 and gens eq {7,29}) or (D eq 6 and N eq 17))
           select 4 else 2;
end function;

curves := GetHyperellipticCandidates();
printf "Loaded %o candidate curves. Verifying %o row(s).\n", #curves, #SEL;

// Group selected rows by star curve (D,N), preserving first-seen order.
group_order := [];                  // ordered list of <D,N> keys
group_rows  := AssociativeArray();  // <D,N> -> [ row indices ]
for i in SEL do
    key := <ROWS[i][1], ROWS[i][2]>;
    if not IsDefined(group_rows, key) then
        group_rows[key] := []; Append(~group_order, key);
    end if;
    Append(~group_rows[key], i);
end for;
printf "Grouped into %o star curve(s).\n\n", #group_order;

// Per-row computed model + status, filled in group by group.
computed := AssociativeArray();   // i -> true/false
curve_of := AssociativeArray();   // i -> CrvHyp (only when computed)
reason_of := AssociativeArray();  // i -> string

// Verdicts accumulate as each group is verified inline (verify-after-each-row),
// so earlier rows report even if a later group hangs or errors.
n_match := 0;
verdicts := [];   // < idx, desc, tag >
models := [* *];  // < idx, desc, C > for every computed model (used in the final biellipticity pass)

// ---- per group: build the models, then verify that group's rows immediately ----
for key in group_order do
    D := key[1]; N := key[2]; idxs := group_rows[key];
    printf "================ star curve X(%o,%o)*  (rows %o) ================\n",
           D, N, [i : i in idxs];
    t0 := Realtime();

    if not exists(Xstar){X : X in curves | X`D eq D and X`N eq N and IsStarCurve(X)} then
        printf "  no star curve found for (D,N)=(%o,%o)\n", D, N;
        for i in idxs do
            computed[i] := false; reason_of[i] := "no star curve found";
            Append(~verdicts, <i, row_desc(i), "NO-MODEL">);
        end for;
        printf "\n"; continue;
    end if;

    imm_W := { curves[j]`W : j in Xstar`CoveredBy };
    imm_idxs  := [i : i in idxs | AllALsFromGens(ROWS[i][3], D*N) in imm_W];
    deep_idxs := [i : i in idxs | AllALsFromGens(ROWS[i][3], D*N) notin imm_W];

    // Immediate covers: one target-restricted EquationsOfCovers for all of them.
    if #imm_idxs gt 0 then
        try
            tgts := { AllALsFromGens(ROWS[i][3], D*N) : i in imm_idxs };
            crv_list, ws, keys := EquationsOfCovers(Xstar, curves : Targets := tgts);
            for i in imm_idxs do
                W := AllALsFromGens(ROWS[i][3], D*N);
                if exists(k){k : k in keys | curves[k]`W eq W} then
                    computed[i] := true; curve_of[i] := crv_list[Index(keys, k)];
                    reason_of[i] := "immediate";
                else
                    computed[i] := false; reason_of[i] := "immediate cover not among computed keys";
                end if;
            end for;
        catch e
            for i in imm_idxs do computed[i] := false; reason_of[i] := Sprintf("ERROR: %o", e`Object); end for;
        end try;
    end if;

    // Deep covers: one build-upward computation for all of them.
    if #deep_idxs gt 0 then
        try
            all_eqns, all_ws := AllEquationsAboveCovers(D, N, curves);
            for i in deep_idxs do
                W := AllALsFromGens(ROWS[i][3], D*N);
                if exists(k){k : k in Keys(all_eqns) |
                             curves[k]`D eq D and curves[k]`N eq N and curves[k]`W eq W} then
                    if #Keys(all_eqns[k]) eq 0 then
                        computed[i] := false; reason_of[i] := "cover not reachable by climbing";
                    else
                        computed[i] := true; reason_of[i] := "climb";
                        curve_of[i] := all_eqns[k][Representative(Keys(all_eqns[k]))];
                    end if;
                else
                    computed[i] := false; reason_of[i] := "no curve id with this W in all_eqns";
                end if;
            end for;
        catch e
            for i in deep_idxs do computed[i] := false; reason_of[i] := Sprintf("ERROR: %o", e`Object); end for;
        end try;
    end if;

    printf "  ---- star curve X(%o,%o)* built in %o s; verifying its row(s) ----\n", D, N, Realtime()-t0;

    // Verify this group's rows NOW (isomorphism to the published table model), so
    // verdicts appear incrementally and a later slow/bad group can't hide them.
    for i in idxs do
        rr := ROWS[i]; h := rr[4]; f := rr[5]; desc := row_desc(i);
        Ctab := HyperellipticCurve(f, h);   // y^2 + h*y = f
        if not computed[i] then
            printf "  [%o] COULD-NOT-COMPUTE: %o\n", desc, reason_of[i];
            Append(~verdicts, <i, desc, "NO-MODEL">);
            continue;
        end if;
        C := curve_of[i];
        Append(~models, <i, desc, C>);
        try
            gC := Genus(C);
            if gC ne 2 then
                printf "  [%o] NOT GENUS 2: computed model has genus %o (%o)\n", desc, gC, reason_of[i];
                Append(~verdicts, <i, desc, Sprintf("NOT-GENUS-2 (g=%o)", gC)>);
            elif IsIsomorphic(C, Ctab) then
                printf "  [%o] MATCH (%o)\n", desc, reason_of[i];
                n_match +:= 1;
                Append(~verdicts, <i, desc, "MATCH">);
            else
                same := G2Invariants(C) eq G2Invariants(Ctab);
                printf "  [%o] MISMATCH (%o); G2Invariants %o\n", desc, reason_of[i],
                       same select "agree (Qbar-isomorphic)" else "differ";
                printf "      computed: y^2 + (%o)*y = %o\n", HyperellipticPolynomials(C);
                Append(~verdicts, <i, desc, same select "MISMATCH(geom-iso)" else "MISMATCH">);
            end if;
        catch e
            printf "  [%o] ERROR during check: %o\n", desc, e`Object;
            Append(~verdicts, <i, desc, "ERROR-CHECK">);
        end try;
    end for;
    printf "\n";
end for;

// Sort verdicts back into row order for the summary.
Sort(~verdicts, func< a, b | a[1] - b[1] >);

printf "\n================ SUMMARY ================\n";
for v in verdicts do
    printf "  %-9o %o\n", v[3], v[2];
end for;
printf "----------------------------------------\n";
printf "  %o / %o matched\n", n_match, #SEL;
printf "========================================\n";

// ---- Biellipticity check over ALL computed models (separate final pass) ----
// A genus-2 curve is bielliptic iff it has an involution besides the
// hyperelliptic one, i.e. iff its geometric automorphism group has order > 2.
// Expected bielliptic (|Aut| = 4): X(14,29)/<w7,w29> and both D=6,N=17 rows;
// every other model is non-bielliptic (|Aut| = 2).
printf "\n========= BIELLIPTICITY CHECK (all models) =========\n";
biell_all_ok := true;
for m in models do
    exp_aut := expected_aut(m[1]);
    exp_bi := exp_aut eq 4;
    try
        if Genus(m[3]) ne 2 then
            printf "  SKIP   %o : model not genus 2\n", m[2];
            biell_all_ok := false;
            continue;
        end if;
        naut := #GeometricAutomorphismGroup(m[3]);
        is_bi := naut gt 2;
        ok := naut eq exp_aut;
        biell_all_ok := biell_all_ok and ok;
        printf "  %-6o %o : |Aut| = %o -> %o (expected %o)\n",
               ok select "OK" else "WRONG", m[2], naut,
               is_bi select "BIELLIPTIC" else "not bielliptic",
               exp_bi select "bielliptic" else "not bielliptic";
    catch e
        printf "  ERROR  %o : %o\n", m[2], e`Object;
        biell_all_ok := false;
    end try;
end for;
printf "---------------------------------------------------\n";
assert biell_all_ok;
printf "  all biellipticity results as expected.\n";
printf "===================================================\n";

exit;
