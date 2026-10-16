// How much has to be pinned EXACTLY to turn each numerical verdict into a proof.
//
// None of these quotients has an Atkin-Lehner overgroup of genus 0, so a hyperelliptic
// involution iota there is exceptional and ITS fixed points are not CM points -- there is no
// quaternionic description of them, which rules out recognising the eight branch points of
// X -> X/iota directly.  What IS exact is the other map: X -> X/u for an Atkin-Lehner u, whose
// branch points ARE the fixed points of u, hence CM points of known discriminant.
//
// The criterion this gives, with iota central so that it descends to an involution of X/u:
//   genus 3, X/u = E of genus 1, four branch points B:
//       X is hyperelliptic  <=>  B splits into two pairs with equal sums in E
//       (a degree-2 map E -> P^1 is P -> c - P, and its fibres are the pairs summing to c)
//   genus 4, X/u = C of genus 2, two branch points B:
//       X is hyperelliptic  <=>  B is stable under the hyperelliptic involution of C,
//       which is unique and always exists, so only two points have to be placed.
AttachSpec("/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/ShimuraQuotients.spec");
SetColumns(0);
cases := [
  <39,11,{1,11,39,429}>, <55,7,{1,35,55,77}>, <6,185,{1,2,15,30,111,185,222,370}>,
  <46,31,{1,2,23,31,46,62,713,1426}>, <6,301,{1,3,14,42,43,129,602,1806}>,
  <1326,1,{1,17,78,1326}>, <6,323,{1,6,19,34,51,114,646,969}>,
  <6,379,{1,2,3,6,379,758,1137,2274}>,
  <6,385,{1,5,6,11,14,21,30,55,66,70,105,154,231,330,770,1155}>,
  <714,5,{1,3,10,14,17,30,35,42,51,105,170,238,510,595,714,1785}>,
  <2002,1,{1,2,13,26,77,154,1001,2002}>, <2030,1,{1,10,14,35,58,145,203,2030}>,
  <6,403,{1,2,13,26,93,186,1209,2418}>, <2418,1,{1,2,13,26,93,186,1209,2418}>,
  <3390,1,{1,5,6,30,113,565,678,3390}>
];
printf "%-10o %-3o %-6o %-5o %-6o %-7o %o\n", "base","g","u","g(X/u)","RH r","CM fix", "what must be pinned";
for c in cases do
  D := c[1]; N := c[2]; W := c[3]; M := D*N;
  g := GenusShimuraCurveQuotient(D, N, W);
  halls := [ m : m in Divisors(M) | GCD(m, M div m) eq 1 ];
  found := false;
  for m in halls do
    if m in W then continue; end if;
    WW := { (a*m) div GCD(a,m)^2 : a in W } join W;
    if #WW ne 2*#W then continue; end if;
    gq := GenusShimuraCurveQuotient(D, N, WW);
    r := (2*g - 2) - 2*(2*gq - 2);
    X := CreateShimuraQuot(D, N, W); X`g := g;
    nfix := CountFixedPointsOnQuotient(m, X);
    printf "%-10o %-3o %-6o %-5o %-6o %-7o %o\n",
       Sprintf("%o_%o", D, N), g, m, gq, r, nfix,
       (g eq 3) select "4 CM points on a genus-1 curve, two pairs of equal sum"
                 else "2 CM points on a genus-2 curve, swapped by its hyperelliptic involution";
    found := true;
  end for;
  if not found then
    printf "%-10o %-3o %-6o %-5o %-6o %-7o %o\n", Sprintf("%o_%o", D, N), g,
       "none", "-", "-", "-", "W = W*: NO Atkin-Lehner map to descend along -- out of reach";
  end if;
end for;
printf "BRANCH DONE\n";
quit;
