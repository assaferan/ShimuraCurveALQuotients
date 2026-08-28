// GAP 1 (even D): does g_p = f_p/2 also hold at p = 2?
// At DN = 330 the prime 2 divides DN, so it is RAMIFIED in some supports and
// the EICHLER prime in others -- exactly the handle we need:
//   c_2^Eich  from the support with s = 2:  ct(theta(165,2))/ct(theta(165,1))
//   c_2^ram   from  ct(theta(2,1))/ct(theta(3,1)) * c_3^ram   (the extra Z/2
//             common to both cancels, and c_3^ram is known in closed form)
// Then test  g_2 := (2+1)/2 c_2^Eich - 1   ==   f_2/2 := (2-1) c_2^ram / 2.
AttachSpec("ShimuraQuotients.spec");
load "vvdata/weyl-campaign/fastcosets.m";
load "vvdata/weyl-campaign/etared.m";
load "vvdata/weyl-campaign/ctlat.m";
PREC := 120; CC := ComplexField(PREC); QSIGN := -1;
DN := 330; lev := 4*DN;
words := [ VVSTWord(g) : g in fastCosetReps(lev) ];
nw := #words;
cv := [ VVWordMatrix(w)[2][1] : w in words ];

// c_2^Eich  (165 = 3*5*11 ramified, 2 Eichler)
a1 := ctTheta(grossGram(165, 1), words, CC, QSIGN);
a2 := ctTheta(grossGram(165, 2), words, CC, QSIGN);
// c_2^ram / c_3^ram
b2 := ctTheta(grossGram(2, 1), words, CC, QSIGN);
b3 := ctTheta(grossGram(3, 1), words, CC, QSIGN);

printf "DN = %o, %o cosets.  Testing g_2 = f_2/2 at the prime 2.\n\n", DN, nw;
vE := {@ @}; vR := {@ @};
n := 0; bad := 0; worst := RealField(20)!0;
for i in [1..nw] do
    if Abs(a1[i]) lt 10^(-60) or Abs(b3[i]) lt 10^(-60) then continue; end if;
    cE := a2[i]/a1[i];
    c3 := (cv[i] mod 3 eq 0) select CC!1 else CC!(-1)/3;   // closed form, odd p
    cR := (b2[i]/b3[i]) * c3;
    Include(~vE, ComplexField(10)!cE); Include(~vR, ComplexField(10)!cR);
    g2 := (CC!3/2)*cE - 1;
    f2 := (CC!1)*cR;
    n +:= 1;
    d := Abs(g2 - f2/2);
    if d gt worst then worst := RealField(20)!d; end if;
    if d gt 10^(-100)*Maximum(Abs(g2),1) then bad +:= 1; end if;
end for;
printf "  c_2^Eich takes %o distinct values: %o\n", #vE,
    [ vE[j] : j in [1..Minimum(5,#vE)] ];
printf "  c_2^ram  takes %o distinct values: %o\n", #vR,
    [ vR[j] : j in [1..Minimum(5,#vR)] ];
printf "\n  g_2 vs f_2/2 : %o cosets, %o mismatches, worst |g_2 - f_2/2| = %o\n",
    n, bad, worst;
printf "  VERDICT: %o\n", bad eq 0 select
    "g_2 = f_2/2 HOLDS -- the proof extends verbatim to even D"
    else "g_2 =/= f_2/2 -- the prime 2 needs its own treatment";
printf "DONE\n";
quit;
