// Subhyperelliptic cover models for X_0(210,13) -- the quotient by W = {1, 2, 21, 42, 65, 130, 1365, 2730} only
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
//
// Built from the Heegner heights of the branch points of X/W -> E, twice: E = X/<W, w_15> = 2730s1
// (branch points of discriminants -312 and -120 over Q(sqrt 2)) and E = X/<W, w_5> = 2730u1
// (-420 and -1092 over Q(i)), by tools/heegner_model.m; see PROVENANCE.md. The third elliptic
// quotient, by w_3, is 210a2. The model is invariant under x -> -1/x.
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[Integers()|1,2,21,42,65,130,1365,2730]] := [* <3, P![ 4, 20, 97, 20, -138, -20, 97, -20, 4 ], P![]> *];
