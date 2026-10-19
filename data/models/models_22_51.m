// Subhyperelliptic cover models for X_0(22,51) -- the quotient by W = {1, 6, 17, 22, 33, 102, 374, 561} only
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0), or
// <genus, "CRV", equations> for a curve given by equations.
//
// Built from the Heegner heights of the branch points of X/W -> E = X/<W, w_187> = 1122f1 (branch
// points of discriminants -4488 and -264, both over Q(sqrt -2)), by tools/heegner_model.m; see
// PROVENANCE.md.
//
// The curve is hyperelliptic, but its canonical image is a conic without a rational point, so it has
// no model y^2 = f(x) over Q. It is stored as the double cover y^2 = A(s, z) of that conic, in
// coordinates s, z of weight 1 and y of weight 2; the residual Atkin-Lehner involution is x -> -x and
// E is the quotient by it, y^2 = A(s, z) on the s-line.
models := AssociativeArray();
models[[Integers()|1,6,17,22,33,102,374,561]] := [* <3, "CRV", [ Strings() |
    "y^2 - 9242176*s^4 + 2871136*s^3*z - 4331508*s^2*z^2 - 124292*s*z^3 - 4787449*z^4",
    "x^2 + 59*s^2 + 118*z^2" ]> *];
