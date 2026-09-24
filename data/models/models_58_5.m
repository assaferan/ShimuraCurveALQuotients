// Subhyperelliptic cover models for X_0(58,5)* -- Guo-Yang / AllEquationsAboveCovers
// models[Sort(W)] := [* <genus, f, h> *] ; model is y^2 + h*y = f (h usually 0).
//
// PARTIAL SET, and deliberately so.  X_0(58,5)* has SEVEN immediate covers, of genera
// 0,1,2,2,3,3,4.  CM-point demand is max(2g+5) over the RETAINED covers, and cmsupply.m
// measured this base SHORT by 3 at the full demand of 13 (supply 10).  Restricting Targets
// to the covers with g <= 2 drops the demand to 9, which the supply meets -- so these are
// the four cover keys reachable at that cap.  The genus-3 entry below lies ABOVE a targeted
// cover rather than being one of the targets.  The dropped high-genus covers need more CM
// points, not a different method.
//
// Produced with IntegralSolution := true (this base's principal parts are non-integral
// under the default arbitrary solution of sol + Kernel) and with the corrected
// slash-constant tolerance in M0MultiplierExact.  It needed all three: without any one of
// them the run fails.
//
// VERIFIED independently by VerifyModelSet: 48 checks, 0 failures -- genus self-consistency,
// genus against the Shimura-curve genus formula, Weil-polynomial divisibility across 3
// nested cover pairs, and trace-formula point counts at p = 3, 7, 11, 13.  None of those
// touches the Borcherds/Schofer path that produced the models, which matters here because
// reaching them required relaxing a guard.
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
models[[ 1, 10, 29, 290 ]] := [* <0, P![ 4/5, 0, 1/5 ], P![]> *];
models[[ 1, 290 ]] := [* <3, P![ 2233, 7366, 10297, 7994, 3788, 1130, 209, 22, 1 ], P![]> *];
models[[ 1, 2, 145, 290 ]] := [* <1, P![ 0, 64/125, 13/125, 32/125, 16/125 ], P![]> *];
models[[ 1, 5, 58, 290 ]] := [* <2, P![ -15, 31, -7, 5, 1 ], P![ 0, 1, 0, 1 ]> *];
