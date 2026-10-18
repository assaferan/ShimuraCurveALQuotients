// For each target: the full Atkin-Lehner group, W, and for every u outside W the genus of
// X/<W,u> (overgroup formed with the Atkin-Lehner group law); also the W-fixed newforms of level
// exactly D*N with their u-eigenvalues (quaternionic: a_p at p | D, -a_p at p | N), so the
// genus-1 choices can be cross-checked against the 1 + 2 splitting.
AttachSpec("/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/ShimuraQuotients.spec");
SetColumns(0);
targets := [
  <55, 7, {1,35,55,77}>,
  <6, 185, {1,2,15,30,111,185,222,370}>,
  <1326, 1, {1,17,78,1326}>,
  <6, 301, {1,3,14,42,43,129,602,1806}>,
  <6, 323, {1,6,19,34,51,114,646,969}>,
  <6, 385, {1,5,6,11,14,21,30,55,66,70,105,154,231,330,770,1155}>,
  <714, 5, {1,3,10,14,17,30,35,42,51,105,170,238,510,595,714,1785}>,
  <46, 31, {1,2,23,31,46,62,713,1426}>,
  <6, 379, {1,2,3,6,379,758,1137,2274}>,
  <6, 649, {1,2,3,6,11,22,33,59,66,118,177,354,649,1298,1947,3894}>
];
for t in targets do
    D := t[1]; N := t[2]; W := t[3]; DN := D*N;
    halls := [ m : m in Divisors(DN) | GCD(m, DN div m) eq 1 ];
    g := GenusShimuraCurveQuotient(D, N, W);
    printf "\n=== X_0^%o(%o)/W, D*N = %o, #AL = %o, #W = %o, genus %o\n", D, N, DN, #halls, #W, g;
    if #W eq #halls then printf "   W is the full group: no u available\n"; continue; end if;
    // newforms
    NFs := Newforms(CuspForms(Gamma0(DN), 2));
    sel := [];
    for orb in NFs do
        f := orb[1];
        if Level(Parent(f)) ne DN then continue; end if;
        dg := Degree(BaseRing(Parent(f)));
        q := qExpansion(f, Max(PrimeDivisors(DN)) + 2);
        ev := AssociativeArray();
        for m in halls do
            e := 1;
            for p in PrimeDivisors(m) do
                ap := Coefficient(q, p);
                e *:= (D mod p eq 0) select ap else -ap;
            end for;
            ev[m] := e;
        end for;
        if &and[ ev[m] eq 1 : m in W ] then Append(~sel, <dg, ev>); end if;
    end for;
    printf "   W-fixed newforms of level exactly %o: %o (sum of degrees %o)\n", DN, #sel, &+[ s[1] : s in sel ];
    seen := {};
    for u in halls do
        if u in W then continue; end if;
        WU := { AtkinLehnerMul(u, m, DN) : m in W } join W;
        if WU in seen then continue; end if;
        Include(~seen, WU);
        gu := GenusShimuraCurveQuotient(D, N, WU);
        nplus := &+[ Integers() | s[1] : s in sel | s[2][u] eq 1 ];
        printf "   u = w_%-5o  coset %o  genus X/<W,u> = %o   newforms with u = +1: %o%o\n",
               u, Sort(SetToSequence(WU diff W)), gu, nplus, (gu eq 1 and nplus eq 1) select "   <-- usable" else "";
    end for;
end for;
print "UCHOICE DONE";
quit;
