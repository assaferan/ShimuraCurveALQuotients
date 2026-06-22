// ============================================================================
// VERIFICATION_genus0.m
//
// Verifies the proposition (paper, "Proposition: rational genus 0"):
//
//   For (m,d) in {(2,5), (5,7), (10,14)},  X(10,7)/<w_m, w_d>  =  P^1_Q.
//   In particular these genus-0 curves all admit rational points.
//
// Strategy (the lowgenus Table-6 / ratpts_table6.m workflow): construct each
// genus-0 cover of the star curve X(10,7)* with the target-restricted
// EquationsOfCovers, confirm it has genus 0, and exhibit a rational point. A
// smooth genus-0 curve over Q with a rational point is isomorphic to P^1_Q, so
// a rational point proves the claimed isomorphism.
//
// Each W = <w_m, w_d> is an index-2 (immediate) cover of X*, so all three are
// solved in a single EquationsOfCovers call. Uses cached polymake solutions
// (level M = 140), so polymake/libnormaliz is not required.
//
// USAGE:  magma VERIFICATION_genus0.m
//
// RUNTIME: ~14 min, ~7 GB peak memory -- almost entirely the one-time
// Borcherds-form build for level M = 140; the three covers themselves are solved
// in a single, cheap EquationsOfCovers call.
// ============================================================================

AttachSpec("ShimuraQuotients.spec");
SetVerbose("ShimuraQuotients", 1);

D := 10; N := 7;
PAIRS := [ {2,5}, {5,7}, {10,14} ];   // the (m,d) generator sets

curves := GetHyperellipticCandidates();
printf "Loaded %o candidate curves.\n", #curves;

assert exists(Xstar){X : X in curves | X`D eq D and X`N eq N and IsStarCurve(X)};

// Solve all three target covers at once.
Targets := { AllALsFromGens(p, D*N) : p in PAIRS };
crv_list, ws, keys := EquationsOfCovers(Xstar, curves : Targets := Targets);

all_ok := true;
verdicts := [];   // < desc, tag >

for p in PAIRS do
    W := AllALsFromGens(p, D*N);
    desc := Sprintf("X(%o,%o)/<%o>", D, N, p);
    printf "\n==== %o ====\n", desc;

    if not exists(k){k : k in keys | curves[k]`W eq W} then
        printf "  cover not among computed keys\n";
        all_ok := false; Append(~verdicts, <desc, "NO-MODEL">); continue;
    end if;
    C := crv_list[Index(keys, k)];

    g := Genus(C);
    if g ne 0 then
        printf "  genus = %o (expected 0)\n", g;
        all_ok := false; Append(~verdicts, <desc, Sprintf("genus %o (expected 0)", g)>); continue;
    end if;

    f := HyperellipticPolynomials(C);
    if Degree(f) le 1 then
        // Already a rational line: trivially P^1_Q.
        printf "  genus 0, model is P^1 (deg f = %o) => isomorphic to P^1_Q\n", Degree(f);
        Append(~verdicts, <desc, "P^1 (deg f <= 1)">);
        continue;
    end if;

    // Genus-0 conic: a rational point gives the isomorphism to P^1_Q.
    con := Conic(C);
    has, pt := HasRationalPoint(con);
    if has then
        printf "  genus 0 with rational point %o => isomorphic to P^1_Q\n", pt;
        printf "    conic: %o\n", con;
        Append(~verdicts, <desc, Sprintf("rational point %o", pt)>);
    else
        printf "  NO rational point on conic %o\n", con;
        all_ok := false; Append(~verdicts, <desc, "NO rational point">);
    end if;
end for;

printf "\n=========== PROPOSITION (rational genus 0) ===========\n";
for v in verdicts do
    printf "  %-22o : %o\n", v[1], v[2];
end for;
printf "------------------------------------------------------\n";
assert all_ok;
printf "  VERIFIED: all three quotients are P^1_Q (admit rational points).\n";
printf "======================================================\n";

exit;
