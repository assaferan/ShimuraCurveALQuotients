// A check the isogeny class CANNOT supply.  The construction builds the model in a coordinate
// where u = w_3 acts as x -> -x, so the four fixed points of u on X sit over x = 0 and x = infinity.
// On  y^2 = kappa * gnorm(mu x^2)  with gnorm(0) = 1 those are y = +-sqrt(kappa) at x = 0, and
// y/x^4 = +-sqrt(kappa*A^3*B^2*C*mu^4) at infinity.  So the model PREDICTS the two quadratic fields
// over which the fixed points of w_3 are defined.  CM theory predicts them independently.
AttachSpec("/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/ShimuraQuotients.spec");
SetColumns(0);
D := 39; N := 11; W := {1,11,39,429};
A := 1391/7688; B := -32364/11449; Cv := -1177/6728; kappa := 13; mu := -403;
printf "from the MODEL:\n";
printf "   fixed points over x = 0        : Q(sqrt(%o))\n", Squarefree(Numerator(kappa)*Denominator(kappa));
lead := kappa*A^3*B^2*Cv*mu^4;
sf := Squarefree(Numerator(lead)*Denominator(lead));
printf "   fixed points over x = infinity : Q(sqrt(%o))   [leading coefficient %o]\n", sf, lead;
printf "\nfrom CM THEORY, the fixed points of each Atkin-Lehner involution on X_0^%o(%o):\n", D, N;
for m in Sort(SetToSequence(W)) do
    if m eq 1 then continue; end if;
    try
        fp := NumFixedPointsByCMOrder(D, N, m);
        printf "   w_%-4o %o\n", m, fp;
    catch e
        printf "   w_%-4o (not available: %o)\n", m, e`Object;
    end try;
end for;
printf "\nand on the quotient by W:\n";
try
    printf "   %o\n", CountFixedPointsOnQuotient(D, N, W, 3);
catch e printf "   (not available)\n"; end try;
print "CERT3911 DONE";
