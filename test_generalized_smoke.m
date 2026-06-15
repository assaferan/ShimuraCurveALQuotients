// Phase-1 smoke test: compile the new package and sanity-check the fixed-point helper.
SetQuitOnError(true);
AttachSpec("ShimuraQuotients.spec");
Attach("GeneralizedComplicatedFixedPoints.m");
print "=== attached + compiled OK ===";

D := 1; N := 72;   // 8 | 72 (V2 available), v3(72)=2 (V3 available), 4 | 72 (S2)
DN := D*N;
printf "X_0(%o,%o): g_X = %o, omega = %o\n", D, N, GenusShimuraCurve(D, N), Omega(DN);

// S2 fixed points on X (vname S2, Q=1 => bare S2)
S2 := Matrix(Integers(),2,2,[2,1,0,2]);
nuS2 := NumFixedPointsNonALOnX(S2, "S2", 1, D, N);
printf "nu(S2 on X_0(%o,%o)) = %o\n", D, N, nuS2;
assert nuS2 ge 0 and IsEven(nuS2);

// V2 fixed points on X
V2 := get_V2(DN);
nuV2 := NumFixedPointsNonALOnX(V2, "V2", 1, D, N);
printf "nu(V2 on X_0(%o,%o)) = %o\n", D, N, nuV2;
assert nuV2 ge 0 and IsEven(nuV2);

// V3 fixed points on X
V3 := get_V3(DN);
nuV3 := NumFixedPointsNonALOnX(V3, "V3", 1, D, N);
printf "nu(V3 on X_0(%o,%o)) = %o\n", D, N, nuV3;
assert nuV3 ge 0 and IsEven(nuV3);

// Exercise the helper enumerators on a trivial W = {1}.
W := {Integers() | 1};
n2s := GeneralizedN2Candidates(D, N, W);
printf "N2 candidates (W={1}) for (%o,%o): %o\n", D, N, n2s;
vs, names, bads := AvailableNonALInvolutions(D, N, W);
printf "available non-AL involutions: %o\n", names;

print "=== smoke test done ===";
quit;
