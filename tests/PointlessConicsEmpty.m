// tests/PointlessConicsEmpty.m
//
// Regression test for EquationsAbovePointlessConics on an EMPTY equation set.
//
// The intrinsic opened with
//     all_keys := Keys(all_eqns);
//     starcurve := Representative(curves[Maximum(all_keys)]`Covers);
// -- Maximum on an unguarded key set.  When EquationsAboveP1s upstream finds no equations at
// all, all_eqns is empty and this throws
//     Runtime error in 'Maximum': Argument 1 is not non-empty
//
// The crash predates the Targets parameter, but restricting Targets makes it much likelier to
// reach: fewer retained covers means a higher chance nothing is found upstream.  X_0(74,3) hit
// it under a g <= 2 cap after clearing the integrality and CM-supply gates, and the failure was
// unlocalised for a while because the driver truncated the error text.
//
// The test is deliberately CHEAP: the guard returns before it touches `curves`, so an empty
// (untyped) sequence suffices and no curve list has to be built.  It therefore costs only
// AttachSpec, and belongs in the CI matrix rather than being excluded as expensive.

printf "PointlessConicsEmpty.m: empty all_eqns must return, not crash...";

pce_e := AssociativeArray();
pce_w := AssociativeArray();

pce_a, pce_b := EquationsAbovePointlessConics(pce_e, pce_w, []);

error if #Keys(pce_a) ne 0,
    Sprintf("PointlessConicsEmpty: expected 0 keys back, got %o", #Keys(pce_a));
error if #Keys(pce_b) ne 0,
    Sprintf("PointlessConicsEmpty: expected 0 ws back, got %o", #Keys(pce_b));

printf "Done!\n";
