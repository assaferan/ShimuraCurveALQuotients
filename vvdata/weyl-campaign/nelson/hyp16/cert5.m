AttachSpec("/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/ShimuraQuotients.spec");
SetColumns(0);
D := 39; N := 11; DN := D*N; W := {1,11,39,429};
coset := { AtkinLehnerMul(3, m, DN) : m in W };
printf "the coset 3*W under the Atkin-Lehner group law: %o\n", Sort(SetToSequence(coset));
tot := 0;
printf "\nfixed points on X_0^%o(%o) of each element, by CM discriminant:\n", D, N;
for m in Sort(SetToSequence(coset)) do
    n := NumFixedPoints(D, N, m); tot +:= n;
    fp := NumFixedPointsByCMOrder(D, N, m);
    ks := Sort(SetToSequence(Keys(fp)));
    printf "   w_%-5o total %-3o  %o\n", m, n, [<d, fp[d]> : d in ks | fp[d] gt 0];
end for;
printf "\nsum %o, divided by #W = %o gives %o fixed points on X/W\n", tot, #W, tot/#W;
printf "Riemann-Hurwitz for X/W (genus 3) -> E (genus 1), degree 2, demands 4.\n";
printf "\nthe quadratic fields involved, against the model's Q(sqrt(13)) and Q(sqrt(-11)):\n";
for m in Sort(SetToSequence(coset)) do
    fp := NumFixedPointsByCMOrder(D, N, m);
    for d in Sort(SetToSequence(Keys(fp))) do
        if fp[d] gt 0 then printf "   w_%-5o disc %-7o field Q(sqrt(%o))\n", m, d, Squarefree(Abs(d))*Sign(d); end if;
    end for;
end for;
print "CERT5 DONE";
