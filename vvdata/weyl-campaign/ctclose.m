// CLOSING GAP A at DN = 1155 (odd, so 2 does not divide DN and the 2-part is
// the lattice's own Z/2 in every case).
//
//   rho          = c_2^L * prod_{p|DN} c_p^{L,p}     (tensor factorisation)
//   c_p^{L,p}    = c_p^ram                            (ramification is LOCAL)
//   RHS          = c_2^G * prod_{p|DN} c_p^ram        (today's algebra)
//   c_p^ram(g)   = 1 if p|c, -1/p otherwise           (today's closed form)
//
// so the identity holds iff c_2^L = c_2^G, and BOTH are now computable by
// dividing out the known odd part.  This is the last link.
AttachSpec("ShimuraQuotients.spec");
load "vvdata/weyl-campaign/fastcosets.m";
load "vvdata/weyl-campaign/etared.m";
load "vvdata/weyl-campaign/ctlat.m";
PREC := 120; CC := ComplexField(PREC); QSIGN := -1;

DN := 1155; lev := 4*DN;
words := [ VVSTWord(g) : g in fastCosetReps(lev) ];
nw := #words;
cvals := [ VVWordMatrix(w)[2][1] : w in words ];
pr := PrimeDivisors(DN);

// the known odd part, in closed form
oddpart := function(i, primes)
    z := CC!1;
    for p in primes do
        z *:= (cvals[i] mod p eq 0) select CC!1 else CC!(-1)/p;
    end for;
    return z;
end function;

// rho for the SHIMURA CURVE lattice, via the same closed form
Ld := ShimuraCurveLattice(DN, 1);
rho := ctTheta(ChangeRing(Ld`Q, Rationals()), words, CC, QSIGN);

// a Gross lattice of the support: theta(D',1), ramified at the primes of D'
s := 3; Dp := DN div s;                 // D' = 385, s = 3
g1 := ctTheta(grossGram(Dp, 1), words, CC, QSIGN);

printf "DN = %o : comparing the 2-adic factors\n", DN;
printf "  c_2^L = rho / prod_{p|DN} c_p^ram        (all four primes ramified in L)\n";
printf "  c_2^G = ct(theta(%o,1)) / prod_{p|%o} c_p^ram\n\n", Dp, Dp;
bad := 0; n := 0; worst := RealField(20)!0;
prD := PrimeDivisors(Dp);
for i in [1..nw] do
    oL := oddpart(i, pr);  oG := oddpart(i, prD);
    if Abs(oL) lt 10^(-60) or Abs(oG) lt 10^(-60) then continue; end if;
    c2L := rho[i]/oL;  c2G := g1[i]/oG;
    n +:= 1;
    d := Abs(c2L - c2G);
    if d gt worst then worst := RealField(20)!d; end if;
    if d gt 10^(-100) then bad +:= 1; end if;
end for;
printf "  %o cosets compared, %o mismatches, worst |c_2^L - c_2^G| = %o\n",
    n, bad, worst;
printf "  VERDICT: %o\n", bad eq 0 select
    "c_2^L = c_2^G  ==>  THE IDENTITY IS PROVEN at this base"
    else "MISMATCH -- the reduction or a local form is wrong";
// and exhibit c_2 itself
vals := {@ @};
for i in [1..Minimum(nw, 4000)] do
    oL := oddpart(i, pr);
    if Abs(oL) gt 10^(-60) then Include(~vals, ComplexField(10)!(rho[i]/oL)); end if;
end for;
printf "  c_2 takes %o distinct values; first few: %o\n", #vals,
    [ vals[j] : j in [1..Minimum(6,#vals)] ];
printf "DONE\n";
quit;
