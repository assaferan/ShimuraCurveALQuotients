// Subhyperelliptic cover models for X_0(6,323) -- the quotient by W = {1, 6, 19, 34, 51, 114, 646, 969} only
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0), or
// <genus, "CRV", equations> for a curve given by equations.
//
// Built from the Heegner heights of the branch points of X/W -> E = X/<W, w_2> = 1938h1, by the
// search described in PROVENANCE.md. The branch locus is a single place of degree 4, of discriminant
// -7752, whose field is the quartic field of discriminant 62016 with quadratic subfield Q(sqrt -2).
//
// The curve is hyperelliptic, but its canonical image is a conic without a rational point, so it has
// no model y^2 = f(x) over Q. It is stored as the double cover y^2 = A(s, z) of that conic, in
// coordinates s, z of weight 1 and y of weight 2; the quotient by the residual Atkin-Lehner
// involution is y^2 = A(s, z) on the s-line, whose Jacobian is 1938h1.
models := AssociativeArray();
models[[Integers()|1,6,19,34,51,114,646,969]] := [* <3, "CRV", [ Strings() |
    "y^2 - 46521*s^4 - 32184*s^3*z - 101968*s^2*z^2 - 91392*s*z^3 - 41984*z^4",
    "x^2 + 22*s^2 + 44*z^2" ]> *];
