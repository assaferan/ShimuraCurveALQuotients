// Step 2 check: numeric eta-quotient evaluation vs the q-expansion.
AttachSpec("ShimuraQuotients.spec");
Attach("vvlib.m");

CC := ComplexField(60);
i := CC.1;

R := EtaQuotientsRing(4, 1);
printf "ds = %o\n", R`ds;
// Delta = eta(z)^24 ; and a mixed one eta(z)^16 eta(2z)^4 eta(4z)^4 (sum d*r = 16+8+16 = 40, not 24|)
tests := [ [24,0,0], [16,4,0], [8,8,2] ];
for exps in tests do
    if (&+[R`ds[k]*exps[k] : k in [1..#exps]]) mod 24 ne 0 then
        printf "exps=%o  SKIP (24 does not divide sum d*r = %o)\n", exps,
               &+[R`ds[k]*exps[k] : k in [1..#exps]];
        continue;
    end if;
    f := EtaQuotient(R, exps);
    ser := qExpansionAtoo(f, 30);
    for y in [0.7, 1.2] do
        z := CC ! (0.13 + y*i);
        q := Exp(2*Pi(CC)*i*z);
        fromser := &+[CC | Coefficient(ser, n)*q^n : n in [Valuation(ser)..25]];
        num := VVEtaEval(f, z);
        printf "exps=%-12o y=%-4o  |numeric - series| = %o\n", exps, y, Abs(num - fromser);
    end for;
end for;
quit;
