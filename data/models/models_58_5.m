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
models[[Integers()|1,10]] := [* <7, P![ -704/625, 4352/625, -4224/125, 61952/625, -132288/625, 196608/625, -24704/125, 20224/625, 243968/625, -20224/625, -24704/125, -196608/625, -132288/625, -61952/625, -4224/125, -4352/625, -704/625 ], P![]> *];
models[[Integers()|1,10,58,145]] := [* <3, P![ -704/625, 9984/625, -76928/625, 383744/625, -1338048/625, 3310592/625, -5630976/625, 6111232/625, -3080192/625 ], P![]> *];
models[[Integers()|1,2,5,10]] := [* <4, P![ -11264/625, 182272/625, -1606656/625, 1880064/125, -39842816/625, 126486528/625, -12123136/25, 542818304/625, -695320576/625, 587464704/625, -49283072/125 ], P![]> *];
models[[Integers()|1,29]] := [* <5, P![ -2816/625, 4608/125, -17664/125, 50688/125, -3328/5, 17408/25, -230912/625, -17408/25, -3328/5, -50688/125, -17664/125, -4608/125, -2816/625 ], P![]> *];
models[[Integers()|1,2,29,58]] := [* <3, P![ -45056/625, 729088/625, -5705728/625, 28819456/625, -100323328/625, 248782848/625, -431845376/625, 97583104/125, -12320768/25 ], P![]> *];
models[[Integers()|1,5,29,145]] := [* <2, P![ -2816/625, 39936/625, -262656/625, 1076224/625, -560896/125, 4558848/625, -770048/125 ], P![]> *];
