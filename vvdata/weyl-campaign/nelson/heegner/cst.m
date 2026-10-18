// Test of the explicit Gross-Zagier constant: r_pred = 2|W| * 2^(mu-2) * u^2 * Area(E)/(deg(pi_opt) * Area(E_opt)),
// mu = number of primes dividing both the conductor and d, on the cases where E has conductor D*N.
SetColumns(0);
load "rows.m";
area := func<E | Abs(Imaginary(ComplexConjugate(Periods(E)[1])*Periods(E)[2]))>;
printf "case        E       d   h(d) fund |W|  mu  u2  degpi  A/Aopt   r_obs     r_pred   obs/pred\n";
for r in rows do
    D, N, nW, lab, d, robs := Explode(r);
    E := EllipticCurve(lab);
    if Conductor(E) ne D*N then continue; end if;
    D0 := FundamentalDiscriminant(d); fund := D0 eq d;
    h := ClassNumber(QuadraticField(d)); hd := ClassNumber(BinaryQuadraticForms(d));
    Eopt := EllipticCurve(Substring(lab, 1, #lab-1) cat "1");
    degpi := ModularDegree(Eopt);
    mu := #[p : p in PrimeDivisors(Conductor(E)) | D0 mod p eq 0];
    u2 := D0 eq -3 select 9 else (D0 eq -4 select 4 else 1);
    if not fund then u2 := 1; end if;
    Ar := area(E)/area(Eopt);
    pred := 2*nW * 2^(mu-2) * u2 * Ar / degpi;
    printf "%3o_%-3o %8o %6o %3o %3o  %2o   %o  %2o  %5o  %6.3o  %9o  %9.5o  %8.4o\n", D, N, lab, d, hd, fund, nW, mu, u2, degpi, Ar, robs, pred, robs/pred;
end for;
