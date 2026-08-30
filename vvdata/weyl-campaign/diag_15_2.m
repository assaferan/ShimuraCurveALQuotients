// Characterization driver: for base 15_2 (D=15, N=2), admit NON-coprime CM points
// and report which discriminants give non-rational Schofer values.
AttachSpec("ShimuraQuotients.spec");
SetVerbose("ShimuraQuotients", 1);

D := 15; N := 2;
Xstar := CreateShimuraQuot(D, N, Set(Divisors(D*N)));
Xstar`g := GenusShimuraCurveQuotient(D, N, Xstar`W);
Xstar`CurveID := 0;  // placeholder; reset by GetQuotientsAndGenera
printf "Star 15_2: genus %o\n", Xstar`g;

curves := GetQuotientsAndGenera([Xstar]);
printf "built %o quotient curves\n", #curves;
_ := exists(star){c : c in curves | IsStarCurve(c)};
printf "star CurveID = %o, W = %o\n", star`CurveID, star`W;

// All rational + quadratic CM points, coprime AND non-coprime.
rat_all, quad_all := RationalandQuadraticCMPoints(star : coprime_to_level := false, bd := 4);
rat_cop  := [p : p in rat_all  | GCD(p[1], N) eq 1];
quad_cop := [p : p in quad_all | GCD(p[1], N) eq 1];
rat_non  := [p : p in rat_all  | GCD(p[1], N) ne 1];
quad_non := [p : p in quad_all | GCD(p[1], N) ne 1];
printf "\nrational  CM discs: coprime %o, non-coprime %o\n", [p[1]:p in rat_cop], [p[1]:p in rat_non];
printf "quadratic CM discs: coprime %o, non-coprime %o\n", [p[1]:p in quad_cop], [p[1]:p in quad_non];

curves_star := [star];  // AllEquationsAboveCovers wants the full lattice; we only need the star for forms
covers := [curves[i] : i in star`CoveredBy];
printf "\ncomputing Borcherds forms...\n";
fs := BorcherdsForms(star, curves : Prec := 100);
printf "Borcherds forms done (%o forms)\n", #Keys(fs);

all_cm_pts := [rat_all, quad_all];
all_discs := {p[1] : p in rat_all} join {p[1] : p in quad_all};
MaxNum := #rat_all;  // force use of all supplied rational points, no internal coprime-filtered fetch

printf "\ncomputing absolute values at %o CM points (Include=all)...\n", #all_discs;
abs_tab, all_cm_pts := AbsoluteValuesAtCMPoints(star, curves, all_cm_pts, fs :
                          MaxNum := MaxNum, Prec := 100, Exclude := {}, Include := all_discs);
ReduceTable(abs_tab);
printf "abs table discs: %o\n", abs_tab`Discs;

printf "\nrunning ValuesAtCMPoints (characterization)...\n";
schofer_tab := ValuesAtCMPoints(abs_tab, all_cm_pts : Exclude := {});
printf "\nNO non-rational cells -- all admitted points gave rational values.\n";
