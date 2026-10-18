// Full prediction from Cai-Shu-Tian Theorem 1.5 (special case 2, c = 1, F = Q) for a branch point P = f(CM point)
// on E, f : X_0^D(N) -> X/W -> E of degree 2|W|:
//   h_abs(P) = (1/2) (1/h^2) sum_{chi admissible} u^2 sqrt|d| 2|W| 2^{#Sigma_D} L'(1,E,chi) / (2 deg(pi_opt) Area(E_opt)),
// Sigma_D = primes dividing both the conductor and d; chi runs over the characters of Pic(O_d) (genus characters
// (d1,d2) here, L(s,E,chi) = L(E^d1,s) L(E^d2,s)); chi is admissible when, at every prime v in Sigma_D, the
// Atkin-Lehner eigenvalue of w_v on E (a_v for v | D, -a_v for v | N) equals chi at the ramified prime above v
// (chi_{d1}(v) if v does not divide d1, else chi_{d2}(v)); otherwise the component is torsion.
SetColumns(0); SetDefaultRealFieldPrecision(20);
load "rows.m";
area := func<E | Abs(Imaginary(ComplexConjugate(Periods(E)[1])*Periods(E)[2]))>;
function Lstar(E)
    L := LSeries(E); v := Evaluate(L, 1);
    if Abs(v) gt 1e-8 then return v, 0; end if;
    v1 := Evaluate(L, 1 : Derivative := 1);
    if Abs(v1) gt 1e-8 then return v1, 1; end if;
    return 0, 2;
end function;
printf "case       E        d  h  |W| SigD  obs/pred  admissible characters (d1,d2) and their eigen-signs\n";
nok := 0; nbad := 0;
for r in rows do
    D, N, nW, lab, d, robs := Explode(r);
    E := EllipticCurve(lab);
    if Conductor(E) ne D*N or FundamentalDiscriminant(d) ne d then continue; end if;
    Cl := ClassGroup(QuadraticField(d)); h := #Cl;
    if Exponent(Cl) gt 2 then continue; end if;
    Eopt := EllipticCurve(Substring(lab, 1, #lab-1) cat "1");
    SigD := [p : p in PrimeDivisors(Conductor(E)) | d mod p eq 0];
    eps := AssociativeArray();
    for v in SigD do eps[v] := (D mod v eq 0) select TraceOfFrobenius(E, v) else -TraceOfFrobenius(E, v); end for;
    facs := [<1, d>];
    for d1 in Divisors(Abs(d)) do for s in [1, -1] do
        e1 := s*d1; e2 := d div e1;
        if e1 lt e2 and e1 ne 1 and IsFundamentalDiscriminant(e1) and IsFundamentalDiscriminant(e2) then Append(~facs, <e1, e2>); end if;
    end for; end for;
    assert #facs eq h;
    u2 := d eq -3 select 9 else (d eq -4 select 4 else 1);
    tot := 0; used := "";
    for f in facs do
        d1 := f[1]; d2 := f[2];
        adm := true; signs := "";
        for v in SigD do
            chiv := (d1 mod v ne 0) select KroneckerSymbol(d1, v) else KroneckerSymbol(d2, v);
            signs cat:= Sprintf("%o%o", v, chiv eq eps[v] select "+" else "-");
            if chiv ne eps[v] then adm := false; end if;
        end for;
        if not adm then continue; end if;
        E1 := d1 eq 1 select E else QuadraticTwist(E, d1); E2 := QuadraticTwist(E, d2);
        v1, o1 := Lstar(E1); v2, o2 := Lstar(E2);
        term := (o1 + o2 eq 1) select v1*v2 else 0;
        tot +:= term; used cat:= Sprintf("(%o,%o)[%o] ", d1, d2, signs);
    end for;
    hpred := (1/2) * (1/h^2) * u2 * Sqrt(Abs(d)) * 2*nW * 2^#SigD * tot / (2 * ModularDegree(Eopt) * area(Eopt));
    base := Evaluate(LSeries(E), 1) * Evaluate(LSeries(QuadraticTwist(E, d)), 1 : Derivative := 1) * Sqrt(Abs(d)) / area(E);
    pred := hpred / base;
    q := pred eq 0 select -1 else robs/pred;
    if Abs(q - 1) lt 1e-6 then nok +:= 1; else nbad +:= 1; end if;
    printf "%3o_%-3o %7o %6o %2o  %o  %o   %8.4o  %o\n", D, N, lab, d, h, nW, #SigD, q, used;
end for;
printf "agree: %o   disagree: %o\n", nok, nbad;
