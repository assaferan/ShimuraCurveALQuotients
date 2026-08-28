// PREDICTIVE TEST of theorem-Eeis.md on bases the campaign has NEVER computed.
//
// The theorem says
//     rho = [ mass_s ct(theta(D',s)) - mass_1 ct(theta(D',1)) ] / (mass_s - mass_1)
// for EVERY D' with omega(D)-1 primes, s = D/D'.  Both sides are computed here
// from the LATTICES ALONE -- no eta pool, no E-table, no least squares, no
// fitting of any kind.  If the theorem is right these agree; nothing is tuned.
AttachSpec("ShimuraQuotients.spec");
load "vvdata/weyl-campaign/fastcosets.m";
load "vvdata/weyl-campaign/etared.m";
load "vvdata/weyl-campaign/ctlat.m";
PREC := 120; CC := ComplexField(PREC); QSIGN := -1;

mass := function(Dps, Rs)
    num := &*[ Integers() | p-1 : p in Dps ];
    m := Rationals()!num / (48 * 2^(#Dps - 1));
    for q in Rs do m *:= Rationals()!(q+1)/2; end for;
    return m;
end function;

test := procedure(D)
    lev := 4*D;
    pr := PrimeDivisors(D);
    printf "\n===== FRESH BASE  D = %o = %o,  N = 1,  level %o =====\n",
        D, pr, lev;
    words := [ VVSTWord(g) : g in fastCosetReps(lev) ];
    nw := #words;
    Ld := ShimuraCurveLattice(D, 1);
    rho := ctTheta(ChangeRing(Ld`Q, Rationals()), words, CC, QSIGN);
    rn := Maximum([ Abs(r) : r in rho ]);
    printf "  %o cosets;  rho computed from the SHIMURA lattice, |rho|_max = %o\n",
        nw, RealField(8)!rn;
    for s in pr do
        Dps := [ p : p in pr | p ne s ];
        Dp := &*Dps;
        m1 := mass(Dps, []); ms := mass(Dps, [s]);
        v1 := ctTheta(grossGram(Dp, 1), words, CC, QSIGN);
        vs := ctTheta(grossGram(Dp, s), words, CC, QSIGN);
        den := ms - m1;
        dev := Maximum([ Abs((CC!ms*vs[i] - CC!m1*v1[i])/CC!den - rho[i]) : i in [1..nw] ]);
        printf "  support (%o,1)+(%o,%o): m1=%-7o ms=%-7o ms-m1=%-6o  max|predicted - rho| = %o  (rel %o)  %o\n",
            Dp, Dp, s, m1, ms, den, RealField(8)!dev, RealField(8)!(dev/rn),
            dev lt 10^(-90)*rn select "PREDICTED" else "*** FAILS ***";
    end for;
end procedure;

for D in [39, 51, 55, 46] do test(D); end for;
test(390);                       // a fresh FOUR-prime base: four supports at once
printf "\nDONE\n";
quit;
