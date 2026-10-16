// The genus-4 exact criterion, and the data it needs.
//
// X = X_0^D(N)/W of genus 4, u the Atkin-Lehner involution with C = X/u of genus 2.  Here W has
// index 2 in the full group, so C is the STAR curve.  If X is hyperelliptic its involution iota
// is unique, hence central and defined over Q, so it commutes with u and descends to an
// involution sigma of C; and X/<iota,u> = (X/iota)/u is rational, so C/sigma has genus 0, which
// at genus 2 forces sigma to be THE hyperelliptic involution of C.  The branch divisor B of
// X -> C is the image of Fix(u), which is iota-stable, so:
//
//     X hyperelliptic  ==>  B is stable under the hyperelliptic involution of C
//                      <=>  [B] - K_C lies in Jac(C)[2]
//                      <=>  the two branch points are exchanged by it (B ~ K_C), or both are
//                           Weierstrass points of C.
//
// That is a NECESSARY condition, exact, and it needs only two points placed.  This script
// produces the two points' CM discriminants -- which is the part the repository can already do.
AttachSpec("/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/ShimuraQuotients.spec");
SetColumns(0);
cases := [
  <2002,1,{1,2,13,26,77,154,1001,2002}>, <2030,1,{1,10,14,35,58,145,203,2030}>,
  <6,403,{1,2,13,26,93,186,1209,2418}>,  <2418,1,{1,2,13,26,93,186,1209,2418}>,
  <3390,1,{1,5,6,30,113,565,678,3390}>
];
for c in cases do
    D := c[1]; N := c[2]; W := c[3]; DN := D*N;
    full := { Integers() | d : d in Divisors(DN) | GCD(d, DN div d) eq 1 };
    g := GenusShimuraCurveQuotient(D, N, W);
    gC := GenusShimuraCurveQuotient(D, N, full);
    X := CreateShimuraQuot(D, N, W); X`g := g;
    printf "\n=== X_0^%o(%o)/W  genus %o over the star curve of genus %o\n", D, N, g, gC;
    coset := [ m : m in Sort(SetToSequence(full)) | not (m in W) ];
    printf "    the non-trivial coset of W in W* is %o\n", coset;
    tot := 0; brk := [];
    for m in coset do
        fp := NumFixedPointsByCMOrder(D, N, m);
        for d in Sort(SetToSequence(Keys(fp))) do
            if fp[d] eq 0 then continue; end if;
            Append(~brk, <m, d, fp[d]>);
            tot +:= fp[d];
        end for;
    end for;
    printf "    fixed points on X_0^%o(%o) of the elements of that coset, by <w_m, disc, count>:\n      %o\n",
       D, N, brk;
    printf "    total %o; divided by #W = %o gives %o branch points on X, "
         * "against %o from Riemann-Hurwitz and %o from the repository's count\n",
       tot, #W, tot/#W, (2*g-2) - 2*(2*gC-2), CountFixedPointsOnQuotient(coset[1], X);
    // the fields of definition of those discriminants on the star curve, which is where the two
    // branch points live
    C := CreateShimuraQuot(D, N, full); C`g := gC;
    discs := Sort(SetToSequence({ t[2] : t in brk }));
    for d in discs do
        ok := true; deg := 0;
        try deg := Degree(AbsoluteField(FieldOfDefinitionOfCMPoint(C, d))); catch e ok := false; end try;
        printf "      discriminant %-8o h(d) = %-5o field of definition on the star: %o\n",
           d, ClassNumber(QuadraticForms(d)), ok select Sprint(deg) else "unavailable";
    end for;
end for;
printf "\nG4BRANCH DONE\n";
quit;
