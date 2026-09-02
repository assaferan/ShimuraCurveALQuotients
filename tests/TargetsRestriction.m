// tests/TargetsRestriction.m
//
// Targets defaults to {}, so without this test no CI job ever executes the restriction path.
//
// Targets exists because CM-point demand is num_vals = max(2g+5) over the RETAINED covers: one
// high-genus sibling inflates the requirement for every cover of the base, and dropping it is
// the documented rescue for a CM-SHORT base.  That is how X_0(58,5)* was brought into range
// (demand 13 -> 9 at g <= 2) and produced the four models in data/models/models_58_5.m.
//
// X_0(33,2) is cheap (~19 s per BorcherdsForms call) and has seven immediate covers of genera
// 0,0,1,1,1,1,1.  Restricting Targets to the two genus-0 covers must:
//   * keep strictly fewer forms  (9 -> 4: the two hauptmoduls, keys -1/-2, are always retained
//     on top of the two targets), and
//   * still produce the two hauptmoduls, since the search needs them regardless of Targets.
//
// Both are asserted.  Counting alone would pass if the restriction silently dropped the
// hauptmoduls too, which would break every downstream consumer of fs[-1] / fs[-2].
//
// NOTE this pins BEHAVIOUR (fewer forms, hauptmoduls kept), not the exact count 4, so it does
// not become brittle if the base's cover structure is ever recomputed.

printf "TargetsRestriction.m: Targets must shrink the cover set and keep the hauptmoduls...";

tr_D := 33; tr_N := 2;
tr_X := CreateShimuraQuot(tr_D, tr_N, Set(Divisors(tr_D*tr_N)));
tr_X`g := GenusShimuraCurveQuotient(tr_D, tr_N, tr_X`W);
tr_X`CurveID := 0;
tr_curves := GetQuotientsAndGenera([tr_X]);
_ := exists(tr_star){c : c in tr_curves | IsStarCurve(c)};

tr_gs := [tr_curves[i]`g : i in tr_star`CoveredBy];
tr_min := Minimum(tr_gs);
tr_sub := {tr_curves[i]`W : i in tr_star`CoveredBy | tr_curves[i]`g eq tr_min};

error if #tr_sub eq #tr_star`CoveredBy,
    "TargetsRestriction: the chosen subset is not a proper restriction; test is vacuous";

tr_all := BorcherdsForms(tr_star, tr_curves : Prec := 100);
tr_res := BorcherdsForms(tr_star, tr_curves : Prec := 100, Targets := tr_sub);

error if #Keys(tr_res) ge #Keys(tr_all),
    Sprintf("TargetsRestriction: restricting Targets did not reduce the form count (%o -> %o)",
            #Keys(tr_all), #Keys(tr_res));

// the two hauptmoduls must survive any restriction -- everything downstream reads fs[-1], fs[-2]
error if not (-1 in Keys(tr_res)) or not (-2 in Keys(tr_res)),
    "TargetsRestriction: a hauptmodul (key -1 or -2) was dropped by the restriction";

printf "Done!\n";
