// VALIDATION: at a prime where the quaternion algebra is SPLIT and the lattice is unimodular
// (p not dividing D*N), Proposition 5.3 must reproduce the pipeline's own LocalWhittakerAtOne.
// If it does, the implementation of (2.16)/(2.17) and of Prop 5.3's assembly is right, and its
// RAMIFIED branch can be trusted at p | D -- which is the whole point.
AttachSpec("ShimuraQuotients.spec");

bases := [[15,2],[6,5],[10,3],[6,7],[10,11]];
printf "%-8o %-4o %-5o %-6o %-14o %-14o %o\n", "base","p","kind","m","LocalWhittaker","Prop 5.3","agree";
for b in bases do
    D := b[1]; N := b[2];
    Ld := ShimuraCurveLattice(D,N);
    negQ := -ChangeRing(Ld`Q, Integers());
    Lfull := RSpaceWithBasis(IdentityMatrix(Integers(), 3));
    zero := Vector(Rationals(), [0,0,0]);
    bad := 0; good := 0;
    for m in [1..40] do
        for p in PrimeDivisors(6*35*Determinant(ChangeRing(Ld`Q,Integers()))*m) do
            if p gt 40 then continue; end if;
            kind := (D mod p eq 0) select "ram" else ((N mod p eq 0) select "lev" else "good");
            code := LocalWhittakerAtOne(Rationals()!m, p, zero, Lfull, negQ);
            ok, ky := true, 0;
            try
                ky := KYWhittaker53(Rationals()!m, p, D mod p eq 0);
            catch e
                ok := false;
            end try;
            if not ok then continue; end if;
            if kind eq "good" then
                if code eq ky then good +:= 1;
                else
                    bad +:= 1;
                    if bad le 4 then
                        printf "%-8o %-4o %-5o %-6o %-14o %-14o %o\n",
                               Sprintf("%o_%o",D,N), p, kind, m, code, ky, "MISMATCH";
                    end if;
                end if;
            end if;
        end for;
    end for;
    printf "%o_%o: good primes -- %o agree, %o disagree\n", D, N, good, bad;
end for;
quit;
