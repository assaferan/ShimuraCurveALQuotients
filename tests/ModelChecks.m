// tests/ModelChecks.m
//
// CI check that every committed subhyperelliptic cover model in data/models/models_D_N.m passes
// the independent verification in ModelVerification.m:
//
//   [1] Genus(model) equals the genus recorded beside it,
//   [2] that genus equals X`g from the Shimura-curve genus formula,
//   [3] Weil-polynomial divisibility across nested cover keys W1 subset W2, and
//   [4] the trace-formula point count: #(X/W)(F_p) from ComputePointsViaTrace (Eichler-Selberg on
//       the W-fixed part of the D-new space) equals the model's actual point count.
//
// None of these uses the Borcherds/Schofer CM machinery that produced the models, so they are
// genuine cross-checks.  Check [4] reuses the codebase's own trace-formula point count (the same
// routine behind the Weil/automorphism filters), so it builds no modular-symbol spaces.
//
// Model files are auto-discovered, so new models are covered without editing this file.
//
// NB: this file is deliberately written as TOP-LEVEL statements rather than wrapped in a
// procedure.  run_tests.m executes tests via `eval`, and an `eval` inside a procedure that also
// closes over an outer variable segfaults Magma 2.29; keeping everything at top level avoids it.
// The model sets are read via `eval` (run_tests.m's own idiom) because Magma's `load` requires a
// literal filename and so cannot iterate over a discovered file list.

model_files := Split(Pipe("ls data/models/models_*.m 2>/dev/null", ""), "\n");
model_files := [f : f in model_files | #f gt 0];
error if #model_files eq 0, "ModelChecks: no model files found under data/models/";

mc_total_chk := 0;
mc_total_fail := 0;
mc_bad := [];

for mc_f in model_files do
    mc_parts := Split(mc_f, "/");
    mc_name  := mc_parts[#mc_parts];              // models_<D>_<N>.m
    mc_core  := mc_name[8..#mc_name-2];           // strip "models_" and ".m"
    mc_dn    := Split(mc_core, "_");
    if #mc_dn ne 2 then continue; end if;
    mc_D := StringToInteger(mc_dn[1]);
    mc_N := StringToInteger(mc_dn[2]);
    models := eval (Read(mc_f) cat "\nreturn models;");
    mc_c, mc_fl := VerifyModelSet(models, mc_D, mc_N : Verbose := false);
    printf "  %o_%o: %o checks, %o failures\n", mc_D, mc_N, mc_c, mc_fl;
    mc_total_chk +:= mc_c;
    mc_total_fail +:= mc_fl;
    if mc_fl ne 0 then Append(~mc_bad, Sprintf("%o_%o", mc_D, mc_N)); end if;
end for;

printf "ModelChecks: %o model file(s), %o checks, %o failures\n", #model_files, mc_total_chk, mc_total_fail;
error if mc_total_chk eq 0, "ModelChecks: no checks ran -- verification is vacuous";
error if mc_total_fail ne 0, Sprintf("ModelChecks: FAILURES in %o", mc_bad);
