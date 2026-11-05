// THE BRANCH LOCUS OF  Y = X_0^2002(1)/<w_2, w_13, w_77>  -->  X_0^2002(1)*  (paper, Section 3, part (1)).
//
//   magma -b tools/proof-2002/stepH-branch.m < /dev/null
//
// A point of the star is a branch point of Y -> star exactly when its orbit upstairs consists of
// fixed points of some w_m with m outside H.  Ogg's formula counts those fixed points by CM
// discriminant; Riemann-Hurwitz for a genus-four cover of a genus-two curve allows exactly two
// branch points, so each discriminant that occurs at all contributes a single rational point.

AttachSpec("ShimuraQuotients.spec");
D := 2002; N := 1;
H := {1, 2, 13, 26, 77, 154, 1001, 2002};
G := { m : m in Divisors(D*N) | GCD(m, (D*N) div m) eq 1 };

printf "genus of Y = X/H: %o;  genus of the star X/G: %o\n",
       GenusShimuraCurveQuotient(D, N, H), GenusShimuraCurveQuotient(D, N, G);
printf "Riemann-Hurwitz for Y -> star: 2*4 - 2 = 2*(2*2 - 2) + r, so r = %o branch points\n\n", 6 - 4;

printf "fixed points of w_m on X_0^%o(%o), m outside H, by discriminant:\n", D, N;
total := 0;
for m in Sort(Setseq(G diff H)) do
    e := NumFixedPointsByCMOrder(D, N, m);
    printf "   w_%-4o %2o fixed points   %o\n", m, NumFixedPoints(D, N, m),
           [ <d, e[d]> : d in Sort(Setseq(Keys(e))) ];
    total +:= NumFixedPoints(D, N, m);
end for;
printf "%o fixed points in all, in two discriminants; two branch points, so one orbit of eight per discriminant\n", total;
quit;
