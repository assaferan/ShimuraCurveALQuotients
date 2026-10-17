// Smallest reproducing case for the GenusShimuraCurveQuotient failure.
// The intrinsic asserts its own result is an integer (ShimuraQuotients.m line 254) and the
// assertion fires, so the genus of the BASE Shimura curve cannot be computed at these
// discriminants -- W = {1}, nothing to do with quotients.
AttachSpec("/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/ShimuraQuotients.spec");
SetColumns(0);
bad := [];
for D in [6..500] do
    if not IsSquarefree(D) then continue; end if;
    if IsOdd(#PrimeDivisors(D)) then continue; end if;   // must be a quaternion discriminant
    ok := true;
    try g := GenusShimuraCurveQuotient(D, 1, {Integers()|1}); catch e ok := false; end try;
    if not ok then Append(~bad, D); end if;
end for;
printf "discriminants D <= 500 (squarefree, even number of prime factors, N = 1, W = {1})\n";
printf "where GenusShimuraCurveQuotient fails its own integrality assertion:\n   %o\n", bad;
if #bad gt 0 then
    D := bad[1];
    printf "\nSMALLEST CASE: D = %o = %o\n", D, Factorisation(D);
    printf "   reproduce with:  GenusShimuraCurveQuotient(%o, 1, {Integers()|1});\n", D;
    // what the pieces are, to localise it
    printf "   arithmetic volume of X_0^%o(1): %o\n", D,
        ArithmeticVolume(FuchsianGroup(MaximalOrder(QuaternionAlgebra(D))));
    for q in [2,3,4,6] do
        try
            e := NumberOfEllipticPointsByCMOrder(CreateShimuraQuot(D,1,{Integers()|1}));
            printf "   elliptic point data available: %o\n", Keys(e);
            break;
        catch er printf "   elliptic point data unavailable for order %o\n", q;
        end try;
    end for;
end if;
printf "G410SCAN DONE\n";
quit;
