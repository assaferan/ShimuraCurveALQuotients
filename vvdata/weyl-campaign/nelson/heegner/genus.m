// Class-number-2 and -4 fundamental cases: r_obs/r_pred should equal (1/h^2) * sum over genus characters
// chi = (d1,d2) of L^(')(E^d1,1) L^(')(E^d2,1) / (L(E,1) L'(E^d,1)), the derivative taken on the factor
// whose value vanishes (odd analytic rank).
SetColumns(0); SetDefaultRealFieldPrecision(20);
load "rows.m";
area := func<E | Abs(Imaginary(ComplexConjugate(Periods(E)[1])*Periods(E)[2]))>;
function Lstar(E)   // L(E,1) if nonzero, else L'(E,1); returns value and order used
    L := LSeries(E); v := Evaluate(L, 1);
    if Abs(v) gt 1e-8 then return v, 0; end if;
    v1 := Evaluate(L, 1 : Derivative := 1);
    if Abs(v1) gt 1e-8 then return v1, 1; end if;
    return 0, 2;
end function;
printf "case       E        d  h   obs/pred   (1/h^2) sum_chi    characters used\n";
for r in rows do
    D, N, nW, lab, d, robs := Explode(r);
    E := EllipticCurve(lab);
    if Conductor(E) ne D*N then continue; end if;
    if FundamentalDiscriminant(d) ne d then continue; end if;
    h := ClassNumber(BinaryQuadraticForms(d));
    if h eq 1 then continue; end if;
    Eopt := EllipticCurve(Substring(lab, 1, #lab-1) cat "1");
    mu := #[p : p in PrimeDivisors(Conductor(E)) | d mod p eq 0];
    pred := 2*nW * 2^(mu-2) * area(E) / (ModularDegree(Eopt) * area(Eopt));
    // genus characters: factorisations d = d1*d2 into fundamental discriminants, d1 < d2 unordered, incl. (1, d)
    facs := [<1, d>];
    for d1 in Divisors(Abs(d)) do
        for s in [1, -1] do
            e1 := s*d1; e2 := d div e1;
            if e1 lt e2 and e1 ne 1 and IsFundamentalDiscriminant(e1) and IsFundamentalDiscriminant(e2) then Append(~facs, <e1, e2>); end if;
        end for;
    end for;
    base := Evaluate(LSeries(E), 1) * Evaluate(LSeries(QuadraticTwist(E, d)), 1 : Derivative := 1);
    tot := 0; used := "";
    for f in facs do
        E1 := f[1] eq 1 select E else QuadraticTwist(E, f[1]); E2 := QuadraticTwist(E, f[2]);
        v1, o1 := Lstar(E1); v2, o2 := Lstar(E2);
        term := (o1 + o2 eq 1) select v1*v2 else 0;
        tot +:= term; used cat:= Sprintf("(%o,%o:%o%o) ", f[1], f[2], o1, o2);
    end for;
    printf "%3o_%-3o %7o %6o %2o  %8.4o   %8.4o    %o\n", D, N, lab, d, h, robs/pred, tot/(h^2*base), used;
end for;
