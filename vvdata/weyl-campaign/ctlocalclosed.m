// Verify the CLOSED FORMS for the local Weil (0,0) entries at odd p | DN:
//     c_p^Eich(g) = 1 if p|c,  1/p  if p not| c     [hyperbolic plane, sig 0]
//     c_p^ram (g) = 1 if p|c, -1/p  if p not| c     [anisotropic plane, sig 4]
// whence c_p^ram = eps_p c_p^Eich with eps_p = +-1, and
//     g_p = (p-1)/2 eps_p c_p^Eich,  f_p = (p-1) eps_p c_p^Eich,  g_p/f_p = 1/2.
AttachSpec("ShimuraQuotients.spec");
load "vvdata/weyl-campaign/fastcosets.m";
load "vvdata/weyl-campaign/etared.m";
load "vvdata/weyl-campaign/ctlat.m";
PREC := 120; CC := ComplexField(PREC); QSIGN := -1;
DN := 1155; lev := 4*DN;
words := [ VVSTWord(g) : g in fastCosetReps(lev) ];
nw := #words;
cvals := [ VVWordMatrix(w)[2][1] : w in words ];

printf "=== c_p^Eich: predicted 1 if p|c, 1/p otherwise ===\n";
for p in PrimeDivisors(DN) do
    Dp := DN div p;
    v1 := ctTheta(grossGram(Dp, 1), words, CC, QSIGN);
    vs := ctTheta(grossGram(Dp, p), words, CC, QSIGN);
    bad := 0; n := 0;
    for i in [1..nw] do
        if Abs(v1[i]) lt 10^(-60) then continue; end if;
        n +:= 1;
        pred := (cvals[i] mod p eq 0) select CC!1 else CC!1/p;
        if Abs(vs[i]/v1[i] - pred) gt 10^(-100) then bad +:= 1; end if;
    end for;
    printf "  p=%2o: %o/%o cosets tested, %o mismatches\n", p, n, nw, bad;
end for;

printf "\n=== c_p^ram / c_q^ram: predicted from  1 (p|c) / -1/p (p not| c) ===\n";
pr := PrimeDivisors(DN);
R := AssociativeArray();
for p in pr do R[p] := ctTheta(grossGram(p, 1), words, CC, QSIGN); end for;
for i in [1..#pr-1] do
    p := pr[i]; q := pr[i+1];
    bad := 0; n := 0;
    for j in [1..nw] do
        if Abs(R[q][j]) lt 10^(-60) then continue; end if;
        n +:= 1;
        cp := (cvals[j] mod p eq 0) select CC!1 else CC!(-1)/p;
        cq := (cvals[j] mod q eq 0) select CC!1 else CC!(-1)/q;
        if Abs(R[p][j]/R[q][j] - cp/cq) gt 10^(-100) then bad +:= 1; end if;
    end for;
    printf "  c_%o^ram / c_%o^ram : %o cosets, %o mismatches\n", p, q, n, bad;
end for;

printf "\n=== lambda = g_p/f_p, after dividing out the (p-independent) 2-part ===\n";
printf "   predicted EXACTLY 1/2 times a p-independent factor\n";
for j in [3, 1731, 6915, 12099] do
    vals := [ ];
    for p in pr do
        cE := (cvals[j] mod p eq 0) select CC!1 else CC!1/p;
        eps := (cvals[j] mod p eq 0) select CC!1 else CC!(-1);
        gp := (CC!((p+1)/2))*cE - 1;
        fp := (CC!(p-1))*eps*cE;
        Append(~vals, gp/fp);
    end for;
    printf "   coset %5o: g_p/f_p over p = %o\n", j, [ ComplexField(10)!v : v in vals ];
end for;
printf "DONE\n";
quit;
