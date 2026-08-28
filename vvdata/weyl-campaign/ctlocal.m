// Extract the LOCAL Weil (0,0) entries and test the local identity
//        g_p = lambda * f_p ,   f_p = (p-1) c_p^ram ,  g_p = (p+1)/2 c_p^Eich - 1
// which is EXACTLY equivalent to D'-independence of E_eis's support.
//
// Both local entries are exact RATIOS, because the Weil rep is multiplicative
// on orthogonal sums and the lattices differ at only one prime:
//    c_p^Eich(g)          = ct(theta(DN/p, p))_g / ct(theta(DN/p, 1))_g
//    c_p^ram/c_q^ram (g)  = ct(theta(p,1))_g     / ct(theta(q,1))_g
AttachSpec("ShimuraQuotients.spec");
load "vvdata/weyl-campaign/fastcosets.m";
load "vvdata/weyl-campaign/etared.m";
load "vvdata/weyl-campaign/ctlat.m";
PREC := 120; CC := ComplexField(PREC); QSIGN := -1;

DN := 1155; lev := 4*DN;
words := [ VVSTWord(g) : g in fastCosetReps(lev) ];
nw := #words;
primes := PrimeDivisors(DN);
printf "DN = %o, primes %o, %o cosets\n\n", DN, primes, nw;

// c_p^Eich, one per prime
cE := AssociativeArray();
for p in primes do
    Dp := DN div p;
    v1 := ctTheta(grossGram(Dp, 1), words, CC, QSIGN);
    vs := ctTheta(grossGram(Dp, p), words, CC, QSIGN);
    cE[p] := [ CC | Abs(v1[i]) gt 10^(-60) select vs[i]/v1[i] else CC!0 : i in [1..nw] ];
end for;
// c_p^ram, up to a common factor: use theta(p,1) directly
cR := AssociativeArray();
for p in primes do cR[p] := ctTheta(grossGram(p, 1), words, CC, QSIGN); end for;

// report the local Eichler entry: is it constant on cusp classes?
printf "c_p^Eich, sampled (should be a simple local quantity):\n";
for p in primes do
    vals := {@ @};
    for i in [1..nw] do
        if Abs(cE[p][i]) gt 10^(-60) then Include(~vals, ComplexField(12)!cE[p][i]); end if;
    end for;
    printf "  p=%2o: %o distinct values; first few: %o\n", p, #vals,
        [ vals[j] : j in [1..Minimum(4, #vals)] ];
end for;

// THE TEST: g_p / [ (p-1) c_p^ram ] must be INDEPENDENT of p.
printf "\nlocal identity test:  g_p / [(p-1) c_p^ram]  for each p, at sample cosets\n";
printf "  (equal across p at each coset  <=>  no re-split rule)\n";
samp := [ i : i in [1..nw] | i mod (nw div 8) eq 3 ];
for i in samp do
    ok := true; vals := [ ];
    for p in primes do
        gp := (CC!((p+1)/2))*cE[p][i] - 1;
        fp := (CC!(p-1))*cR[p][i];
        if Abs(fp) lt 10^(-60) then ok := false; break; end if;
        Append(~vals, gp/fp);
    end for;
    if not ok then continue; end if;
    spread := Maximum([ Abs(vals[j] - vals[1]) : j in [1..#vals] ]);
    scale  := Maximum([ Abs(v) : v in vals ]);
    printf "  coset %5o: lambda = %o   spread = %o  (rel %o)\n", i,
        ComplexField(12)!vals[1], RealField(8)!spread,
        RealField(8)!(scale gt 0 select spread/scale else 0);
end for;
printf "DONE\n";
quit;
