// Is the NUMERATOR of the identity E_eis = (Theta_s - Theta_1)/(m_s - m_1)
// independent of D'?
//
// Theta = mass * theta (the UNNORMALIZED mass-weighted genus sum), so the
// claim is that  m_s*ct(theta(D',s)) - m_1*ct(theta(D',1))  is the SAME
// VECTOR for every three-prime D' with D'*s = DN.
//
// The denominator m_s - m_1 is already PROVEN constant: it equals
// prod_{p|D'}(p-1)*(s-1) / (2*48*2^(w-1)), and prod_{p|D'}(p-1)*(s-1)
// telescopes to prod_{p|DN}(p-1) because D' and s partition the primes of DN.
//
// At the S-cusp the numerator is also constant BY HAND: ct = e(sig/8)/sqrt|A|
// with |A| = 2(D'R)^2, so the difference is
//     -e(sig/8) * prod_{p|DN}(p-1) / (384*sqrt2*DN),
// by the same telescoping.  This checks all the other cusps.
AttachSpec("ShimuraQuotients.spec");
load "vvdata/weyl-campaign/fastcosets.m";
load "vvdata/weyl-campaign/etared.m";
load "vvdata/weyl-campaign/ctlat.m";

PREC := 120; CC := ComplexField(PREC);
QSIGN := -1;

mass := function(Dps, Rs)
    num := &*[ Integers() | p-1 : p in Dps ];
    m := Rationals()!num / (48 * 2^(#Dps - 1));
    for q in Rs do m *:= Rationals()!(q+1)/2; end for;
    return m;
end function;

procedure runbase(DN, prime_sets)
    lev := 4*DN;
    words := [ VVSTWord(g) : g in fastCosetReps(lev) ];
    nw := #words;
    printf "\n===== DN = %o, form level %o, %o cosets =====\n", DN, lev, nw;
    ref := [ CC | ]; refname := "";
    for Dps in prime_sets do
        allp := PrimeDivisors(DN);
        s := [ p : p in allp | p notin Dps ][1];
        Dp := &*Dps;
        m1 := mass(Dps, []); ms := mass(Dps, [s]);
        c1 := ctTheta(grossGram(Dp, 1), words, CC, QSIGN);
        cs := ctTheta(grossGram(Dp, s), words, CC, QSIGN);
        num := [ CC | (CC!ms)*cs[i] - (CC!m1)*c1[i] : i in [1..nw] ];
        if #ref eq 0 then
            ref := num; refname := Sprintf("(%o,%o)", Dp, s);
            printf "  D'=%4o s=%2o  m1=%-8o ms=%-8o ms-m1=%-6o  [reference]\n",
                Dp, s, m1, ms, ms-m1;
        else
            dev := Maximum([ Abs(num[i] - ref[i]) : i in [1..nw] ]);
            nrm := Maximum([ Abs(ref[i]) : i in [1..nw] ]);
            printf "  D'=%4o s=%2o  m1=%-8o ms=%-8o ms-m1=%-6o  max|num - num%o| = %o  (rel %o)\n",
                Dp, s, m1, ms, ms-m1, refname,
                RealField(8)!dev, RealField(8)!(dev/nrm);
        end if;
    end for;
    // and the value at the S-cusp, against the hand formula
    prodp := &*[ Integers() | p-1 : p in PrimeDivisors(DN) ];
    printf "  prod_{p|DN}(p-1) = %o ; predicted S-cusp numerator magnitude = %o\n",
        prodp, RealField(12)!(prodp / (384*Sqrt(RealField(30)!2)*DN));
    printf "  measured at the S-coset (word 2): %o\n", RealField(12)!Abs(ref[2]);
end procedure;

runbase(1155, [ [5,7,11], [3,5,7], [3,7,11], [3,5,11] ]);
runbase(330,  [ [3,5,11], [2,3,5], [2,3,11], [2,5,11] ]);
printf "DONE\n";
quit;
