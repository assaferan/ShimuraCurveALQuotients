// For each of the sixteen numerically-hyperelliptic quotients: the full Atkin-Lehner group,
// the genus of X/W, and the genus of every index-2 overgroup of W inside W*.
// If some overgroup has genus 0 then X/W is hyperelliptic BY AN ATKIN-LEHNER INVOLUTION and the
// equation is a double cover of a rational curve.  If none does, the hyperelliptic involution is
// not Atkin-Lehner, which for these groups would contradict the normaliser being W*.
AttachSpec("/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/ShimuraQuotients.spec");
SetColumns(0);

cases := [
  <39, 11, {1,11,39,429}>,
  <55, 7, {1,35,55,77}>,
  <6, 185, {1,2,15,30,111,185,222,370}>,
  <46, 31, {1,2,23,31,46,62,713,1426}>,
  <6, 301, {1,3,14,42,43,129,602,1806}>,
  <1326, 1, {1,17,78,1326}>,
  <6, 323, {1,6,19,34,51,114,646,969}>,
  <6, 379, {1,2,3,6,379,758,1137,2274}>,
  <6, 385, {1,5,6,11,14,21,30,55,66,70,105,154,231,330,770,1155}>,
  <6, 649, {1,2,3,6,11,22,33,59,66,118,177,354,649,1298,1947,3894}>,
  <714, 5, {1,3,10,14,17,30,35,42,51,105,170,238,510,595,714,1785}>,
  <2002, 1, {1,2,13,26,77,154,1001,2002}>,
  <2030, 1, {1,10,14,35,58,145,203,2030}>,
  <6, 403, {1,2,13,26,93,186,1209,2418}>,
  <2418, 1, {1,2,13,26,93,186,1209,2418}>,
  <3390, 1, {1,5,6,30,113,565,678,3390}>
];

for c in cases do
  D := c[1]; N := c[2]; W := c[3];
  M := D*N;
  // full Atkin-Lehner group: all Hall divisors of D*N
  halls := [ m : m in Divisors(M) | GCD(m, M div m) eq 1 ];
  printf "\n=== X_0^%o(%o)   D*N = %o   #W* = %o\n", D, N, M, #halls;
  if #W eq 0 then
     printf "   (W not supplied here)\n";
     continue;
  end if;
  g := GenusShimuraCurveQuotient(D, N, W);
  printf "   W = %o  genus %o   index of W in W* = %o\n", Sort(SetToSequence(W)), g, #halls div #W;
  // index-2 overgroups of W inside W*: W join m*W for m in W* \ W, keeping those of size 2#W
  seen := {};
  for m in halls do
    if m in W then continue; end if;
    WW := { (a*m) div GCD(a,m)^2 : a in W } join W;
    if #WW ne 2*#W then continue; end if;
    if WW in seen then continue; end if;
    Include(~seen, WW);
    gg := GenusShimuraCurveQuotient(D, N, WW);
    printf "      overgroup by w_%o : %o  genus %o%o\n", m, Sort(SetToSequence(WW)), gg,
           gg eq 0 select "   <== GENUS 0: hyperelliptic by Atkin-Lehner" else "";
  end for;
  gstar := GenusShimuraCurveQuotient(D, N, Set(halls));
  printf "   star genus %o\n", gstar;
end for;
printf "\nSTRUCT DONE\n";
quit;
