AttachSpec("/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/ShimuraQuotients.spec");
SetColumns(0);
D := 39; N := 11;
// the involution acting on X/W as "w_3" is the coset 3*W, W = {1,11,39,429}
for m in [3, 33, 117, 1287] do
    fp := NumFixedPointsByCMOrder(D, N, m);
    ks := Sort(SetToSequence(Keys(fp)));
    printf "w_%-5o fixed points by CM discriminant: %o\n", m, [<d, fp[d]> : d in ks];
end for;
printf "\nquadratic fields of those discriminants, and the two the model predicts:\n";
for d in [-3,-33,-117,-1287,-11,-143,-4,-8] do
    printf "   d = %-6o Q(sqrt(%o))\n", d, Squarefree(Abs(d)) * Sign(d);
end for;
printf "\nmodel predicts Q(sqrt(13)) and Q(sqrt(-11)); note 13*(-11) = %o and -11*-143 = %o = 11^2 * 13\n", 13*(-11), (-11)*(-143);
print "CERT2 DONE";
