// The certification checklist (a)-(d) for a genus-3 model produced by hypfinish3.m, for any
// base.  The curve is rebuilt from the three invariants and the two twist constants, so nothing
// is retyped from a printout:
//     y^2 = kappa * gnorm(mu x^2),   gnorm(t) = 1 + t + A t^2 + A^2 B t^3 + A^3 B^2 C t^4,
// with u = w_UM acting as x -> -x and X/<W,u> : Y^2 = kappa gnorm(t).
//
//   magma -b DISC:=39 LEV:=11 WS:=1,11,39,429 UM:=3 AA:=1391/7688 BB:=-32364/11449 \
//            CC0:=-1177/6728 KAP:=13 MU:=-403 certgen.m < /dev/null
//
// (a) the four fixed points of u: CM theory counts the fixed points of the coset u*W (group law
//     AtkinLehnerMul) by CM discriminant, their sum over #W must be 4, and the field of
//     definition of each pair on X/W must match what the model says: the pair over x = 0 is
//     defined over Q(sqrt(kappa)), the pair over x = infinity over Q(sqrt(leading coefficient)).
// (b) the Jacobian of X/<W,u> (Magma's Jacobian of the quartic) against the W-fixed newform of
//     level exactly D*N with quaternionic u = +1, a_p for every good p up to 100.
// (c) the W-fixed newform space of level exactly D*N splits 1 + 2 under u.
// (d) the full numerator of the zeta function of X/W over F_p for good p up to ZCAP (default
//     60), from the trace-formula counts over F_p, F_p^2, F_p^3; and of X/<W,u> likewise.
AttachSpec("/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/ShimuraQuotients.spec");
SetColumns(0);
D := StringToInteger(DISC); N := StringToInteger(LEV); DN := D*N;
W := { StringToInteger(x) : x in Split(WS, ",") };
um := StringToInteger(UM);
A := StringToRational(AA); B := StringToRational(BB); Cv := StringToRational(CC0);
kap := StringToRational(KAP); mu := StringToRational(MU);
ZCAP := StringToInteger(assigned ZC select ZC else "60");
PQ<x> := PolynomialRing(Rationals());
gnorm := 1 + x + A*x^2 + A^2*B*x^3 + A^3*B^2*Cv*x^4;
foct := kap*Evaluate(gnorm, mu*x^2);
fE := kap*gnorm;
WU := { AtkinLehnerMul(um, m, DN) : m in W } join W;
coset := WU diff W;
printf "=== X_0^%o(%o)/%o, u = w_%o, coset %o\n", D, N, Sort(SetToSequence(W)), um, Sort(SetToSequence(coset));
printf "model (u : x -> -x):  y^2 = %o\n", foct;
X3 := HyperellipticCurve(foct);
printf "genus %o;  reduced minimal model: %o\n", Genus(X3), ReducedMinimalWeierstrassModel(X3);
nfail := 0;

// ---------------- (a) fixed points of u ----------------
printf "\n(a) FIXED POINTS OF u\n";
sfk := Squarefree(Numerator(kap)*Denominator(kap));
lead := LeadingCoefficient(foct);
sfl := Squarefree(Numerator(lead)*Denominator(lead));
printf "   model: pair over x = 0 defined over Q(sqrt(%o)), pair over x = infinity over Q(sqrt(%o))\n", sfk, sfl;
predicted := {* sfk, sfl *};
tot := 0; discs := [];
for m in Sort(SetToSequence(coset)) do
    n := NumFixedPoints(D, N, m); tot +:= n;
    fp := NumFixedPointsByCMOrder(D, N, m);
    ks := [ d : d in Sort(SetToSequence(Keys(fp))) | fp[d] gt 0 ];
    printf "   w_%-6o fixed points on X_0^%o(%o): %-4o by discriminant %o\n", m, D, N, n, [ <d, fp[d]> : d in ks ];
    for d in ks do Append(~discs, <d, fp[d]>); end for;
end for;
printf "   total %o over #W = %o: %o fixed points on X/W (Riemann-Hurwitz needs 4)\n", tot, #W, tot/#W;
if tot ne 4*#W then printf "   ** FAIL: not 4 branch points\n"; nfail +:= 1; end if;
XQ := CreateShimuraQuot(D, N, W); XQ`g := 3;
found := {* *};
for dc in discs do
    d := dc[1];
    deg := DegreeOfFieldOfDefinitionOfCMPoint(XQ, d);
    F := FieldOfDefinitionOfCMPoint(XQ, d);
    npts := dc[2] / #W;   // this discriminant's share of the 4 fixed points, on X/W
    if Type(F) eq FldRat then
        printf "   disc %-7o %o point(s) on X/W, field Q (degree %o)\n", d, npts, deg;
        if npts eq 2 then Include(~found, 1); end if;
    elif Degree(F) eq 2 then
        s := Squarefree(Discriminant(F));
        printf "   disc %-7o %o point(s) on X/W, field Q(sqrt(%o)) (degree %o)\n", d, npts, s, deg;
        if npts eq 2 then Include(~found, s); end if;
    else
        printf "   disc %-7o %o point(s) on X/W, field of degree %o: %o\n", d, npts, Degree(F), F;
    end if;
end for;
printf "   CM theory: %o   model: %o   -> %o\n", found, predicted,
       found eq predicted select "AGREE" else "** DISAGREE";
if found ne predicted then nfail +:= 1; end if;

// ---------------- (b), (c) newforms ----------------
printf "\n(b),(c) NEWFORMS OF LEVEL %o FIXED BY W\n", DN;
ps := [ p : p in PrimesUpTo(100) | DN mod p ne 0 ];
plusdim := 0; minusdim := 0; plusf := [];
for orb in Newforms(CuspForms(Gamma0(DN), 2)) do
    f := orb[1];
    if Level(Parent(f)) ne DN then continue; end if;
    dg := Degree(BaseRing(Parent(f)));
    q := qExpansion(f, Max(PrimeDivisors(DN)) + 2);
    ok := true; eu := 0;
    for m in Sort(SetToSequence(WU)) do
        if m eq 1 then continue; end if;
        e := 1;
        for p in PrimeDivisors(m) do e *:= (D mod p eq 0) select Coefficient(q,p) else -Coefficient(q,p); end for;
        if m in W and e ne 1 then ok := false; end if;
        if m eq um then eu := e; end if;
    end for;
    if not ok then continue; end if;
    printf "   W-fixed orbit of degree %o, quaternionic w_%o eigenvalue %o\n", dg, um, eu;
    if eu eq 1 then plusdim +:= dg; Append(~plusf, f); else minusdim +:= dg; end if;
end for;
printf "   (c) splitting under u: %o + %o   -> %o\n", plusdim, minusdim,
       (plusdim eq 1 and minusdim eq 2) select "1 + 2 as required" else "** NOT 1 + 2";
if not (plusdim eq 1 and minusdim eq 2) then nfail +:= 1; end if;
E := Jacobian(GenusOneModel(fE));
E := MinimalModel(E);
printf "   (b) Jacobian of X/<W,w_%o> = Y^2 = %o :\n       %o, conductor %o", um, fE, E, Conductor(E);
try printf ", Cremona %o", CremonaReference(E); catch e ; end try;
printf "\n";
if Conductor(E) ne DN then printf "   ** FAIL: conductor is not %o\n", DN; nfail +:= 1; end if;
if #plusf eq 1 then
    q := qExpansion(plusf[1], 101);
    af := [ Integers() ! Coefficient(q, p) : p in ps ];
    aE := [ TraceOfFrobenius(E, p) : p in ps ];
    printf "   newform a_p, p in %o:\n       %o\n   Jacobian a_p:\n       %o\n", ps, af, aE;
    printf "   -> %o at %o primes\n", af eq aE select "IDENTICAL" else "** DIFFER", #ps;
    if af ne aE then nfail +:= 1; end if;
end if;

// ---------------- (d) zeta numerators ----------------
printf "\n(d) FULL ZETA NUMERATORS against the trace formula, good p <= %o\n", ZCAP;
EQ := CreateShimuraQuot(D, N, WU); EQ`g := 1;
RT<T> := PolynomialRing(Rationals());
procedure zetacheck(C, Xq, g, name, ~nchk, ~nf)
    for p in [ p : p in PrimesUpTo(ZCAP) | (2*DN) mod p ne 0 ] do
        den := LCM([ Denominator(c) : c in Coefficients(HyperellipticPolynomials(C)) ]);
        if den mod p eq 0 then printf "   %o p=%o: model has p in a denominator, skipped\n", name, p; continue; end if;
        Cp := ChangeRing(C, GF(p));
        if not IsNonsingular(Cp) then printf "   %o p=%o: bad reduction of this model, skipped\n", name, p; continue; end if;
        if g eq 1 then
            // Magma's ZetaFunction on a genus-1 quartic goes through EllipticCurve, which needs
            // a rational point; the Jacobian has the same L-polynomial.
            Lmodel := 1 - TraceOfFrobenius(E, p)*T + p*T^2;
        else
            Lmodel := RT ! Numerator(ZetaFunction(Cp));
        end if;
        cnt := [ ComputePointsViaTrace(Xq, p, r) : r in [1..g] ];
        s := [ p^r + 1 - cnt[r] : r in [1..g] ];
        e := [ Rationals() | ];
        for r in [1..g] do
            acc := s[r];
            for j in [1..r-1] do acc +:= (-1)^j * e[j] * s[r-j]; end for;
            Append(~e, (-1)^(r-1)*acc/r);
        end for;
        a := [ Rationals() | 1 ] cat [ (-1)^r*e[r] : r in [1..g] ];
        // a_{2g-r} = p^(g-r) a_r, appended in the order a_{g+1}, ..., a_{2g}
        for r in [g-1..0 by -1] do Append(~a, p^(g-r)*a[r+1]); end for;
        Ltrace := &+[ a[i+1]*T^i : i in [0..2*g] ];
        nchk +:= 1;
        if Ltrace ne Lmodel then
            nf +:= 1; printf "   %o p=%o: ** MISMATCH  model %o  trace %o\n", name, p, Lmodel, Ltrace;
        else
            printf "   %o p=%o: agrees  %o\n", name, p, Lmodel;
        end if;
    end for;
end procedure;
nchk := 0; nfz := 0;
C3 := HyperellipticCurve(foct); dd := LCM([ Denominator(c) : c in Coefficients(foct) ]);
C3 := HyperellipticCurve(dd^2*foct);
zetacheck(C3, XQ, 3, "X/W   ", ~nchk, ~nfz);
CE := HyperellipticCurve(fE); dd := LCM([ Denominator(c) : c in Coefficients(fE) ]);
CE := HyperellipticCurve(dd^2*fE);
zetacheck(CE, EQ, 1, "X/<W,u>", ~nchk, ~nfz);
printf "   zeta numerators checked: %o, mismatches: %o\n", nchk, nfz;
nfail +:= nfz;
printf "\nCERTGEN %o: %o failures\n", nfail eq 0 select "ALL CHECKS PASS" else "FAILED", nfail;
print "CERTGEN DONE";
quit;
