// A CM point that the caller asks to keep (Keep, and the hauptmoduls' zeros and poles) must be in the
// table even when it has degree two and the rational CM points do not fill the table, so that the search
// goes on to further points.  The table must then have the number of points asked for, each once.
//
// X_0^22(1)* has 12 rational CM points among the candidates; -91 and -115 are points of degree two
// (the field of definition is quadratic: special_fiber_cm.m stores their hauptmodul values as quadratic
// minimal polynomials).  With 16 points asked for, both must be kept.

printf "Testing that kept CM points of degree two stay in the table...";
curves := GetHyperellipticCandidates();
Xstar := rep{X : X in curves | X`D eq 22 and X`N eq 1 and IsStarCurve(X)};
cand := CandidateDiscriminants(Xstar, curves);
assert -91 in {p[1] : p in cand[2]} and -115 in {p[1] : p in cand[2]};
assert #cand[1] lt 16;
tab := ValuesAtCMPoints(Xstar, curves : MaxNum := 16, Prec := 100, Keep := {-91, -115});
assert -91 in tab`Discs and -115 in tab`Discs;
assert #tab`Discs eq 16 and #Set(tab`Discs) eq 16;
printf " ok\n";
