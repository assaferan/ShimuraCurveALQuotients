// tests/ConicClasses.m
//
// THE FIRST CHECK OF ANY KIND ON THE TWIST CLASS OF A GENUS-0 MODEL ENTRY.
//
// WHY THIS EXISTS.  ModelChecks' four tests are STRUCTURALLY BLIND to a quadratic twist at genus 0
// (see the banner in tests/ModelChecks.m): a conic and its non-square twist share the genus, the
// genus formula, the trivial Weil polynomial, AND the point count over every F_p -- every smooth
// conic over a finite field is isotropic, so both twists are P^1 with exactly p+1 points.  Measured
// 2026-09-12: 282 of 822 committed entries are genus 0, across 75 model files, and NOTHING
// validated their twist class.  That is the largest exempted class in the repo, and it is why
// models_10_3.m's [1,2] drift (3 non-isomorphic genus-0 conics, found by ModelRegen) was invisible
// to CI.
//
// WHAT IT CHECKS, and why it needs no theory.  The entries stored under ONE (D, N, W) key are the
// same quotient X/W expressed over DIFFERENT BASE covers, so they must be isomorphic -- hence must
// share one class in Br(Q)[2].  For y^2 = f(x) with deg f = 2, completing the square gives
// y^2 - a u^2 = -disc/(4a), and since 1/(4a) = a modulo squares the class is the quaternion algebra
//      (a, disc),      a = leading coeff,   disc = b^2 - 4ac.
// So the check is: within each key, all genus-0 entries must give the same RamifiedPrimes.  It
// requires no prediction of WHICH conic is correct -- that is the open arbiter question (Ogg's
// real-points criterion is the likely frame) -- only that the file agrees with itself.
//
// ⚠ WHAT IT DOES *NOT* CATCH.  A key whose entries are ALL wrong in the same way passes.  In
// particular it does NOT resolve models_10_3.m: that drift is committed-vs-REGENERATED, and the
// committed entries agree with each other (all split).  This is an internal-consistency check, not
// an oracle.
//
// ⚠ NON-VACUOUS, AND CHECKED TO BE.  Negative control run 2026-09-12: twisting ONE entry of the
// real 10_3 [1,2] key by the non-square -2 moves its class from [] to [5] and the check FIRES.
// A passing check is not evidence until you know it could have failed.
//
// NB: top-level statements, not a procedure -- run_tests.m executes tests via `eval`, and an `eval`
// inside a procedure that closes over an outer variable segfaults Magma 2.29 (same trap as
// ModelChecks.m and ModelRegen.m).  No pipeline is run: this is arithmetic on committed data.

cc_files := Split(Pipe("ls data/models/models_*.m 2>/dev/null", ""), "\n");
cc_files := [f : f in cc_files | #f gt 0];
error if #cc_files eq 0, "ConicClasses: no model files found under data/models/";

cc_P<ccx> := PolynomialRing(Rationals());
cc_conics := 0;      // genus-0 entries that are genuine conics (deg f = 2)
cc_keys   := 0;      // keys carrying MORE THAN ONE such entry -- the ones that carry a claim
cc_bad    := [];

for cc_f in cc_files do
    cc_parts := Split(cc_f, "/");
    cc_name  := cc_parts[#cc_parts];
    cc_core  := cc_name[8..#cc_name-2];
    cc_dn    := Split(cc_core, "_");
    if #cc_dn ne 2 then continue; end if;
    cc_D := StringToInteger(cc_dn[1]);
    cc_N := StringToInteger(cc_dn[2]);
    models := eval (Read(cc_f) cat "\nreturn models;");

    for cc_k in Keys(models) do
        cc_rams := [];
        for cc_e in models[cc_k] do
            if cc_e[1] ne 0 then continue; end if;          // genus 0 only
            cc_fo := cc_e[2];
            if Type(cc_fo) eq MonStgElt then continue; end if;   // CRV entry: no plane model here
            cc_poly := cc_P ! cc_fo;
            // ⚠ AN ENTRY MAY BE <genus, f, h>, MEANING y^2 + h y = f.  Dropping h gives a DIFFERENT
            // curve, so complete the square: (y + h/2)^2 = f + h^2/4.
            if #cc_e ge 3 and Type(cc_e[3]) ne MonStgElt then
                cc_poly := cc_poly + (cc_P ! cc_e[3])^2/4;
            end if;
            if Degree(cc_poly) ne 2 then continue; end if;  // deg < 2 is already P^1; deg > 2 is not genus 0
            cc_a := Coefficient(cc_poly, 2);
            cc_disc := Coefficient(cc_poly, 1)^2 - 4*cc_a*Coefficient(cc_poly, 0);
            if cc_a eq 0 or cc_disc eq 0 then continue; end if;   // degenerate: not a smooth conic
            cc_conics +:= 1;
            Append(~cc_rams, Sort(RamifiedPrimes(QuaternionAlgebra<Rationals() | cc_a, cc_disc>)));
        end for;
        if #cc_rams le 1 then continue; end if;             // a single entry asserts nothing here
        cc_keys +:= 1;
        if #Set(cc_rams) gt 1 then
            Append(~cc_bad, Sprintf("%o_%o W=%o rams=%o", cc_D, cc_N,
                                    Sort([Integers()| w : w in cc_k]), cc_rams));
        end if;
    end for;
end for;

printf "  %o genus-0 conic(s); %o multi-entry key(s) carry a consistency claim\n",
       cc_conics, cc_keys;

// ⚠ COUNT GUARDS, so a future regression cannot quietly reduce what is compared -- the failure mode
// this repo has hit repeatedly is a check that silently stops checking.  Raise these when the model
// set grows; never lower them to make a run pass.
error if cc_conics lt 215,
      Sprintf("ConicClasses: only %o conics examined, expected >= 215 -- did entry parsing regress?",
              cc_conics);
error if cc_keys lt 38,
      Sprintf("ConicClasses: only %o multi-entry keys, expected >= 38 -- the check may be vacuous",
              cc_keys);

error if #cc_bad gt 0,
      Sprintf("ConicClasses: %o key(s) hold genus-0 entries with DIFFERENT Brauer classes, so they "
              * "cannot all be the same quotient:\n  %o", #cc_bad, cc_bad);

printf "  all %o multi-entry key(s) internally consistent\n", cc_keys;
