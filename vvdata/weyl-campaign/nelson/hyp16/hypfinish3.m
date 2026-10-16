// From the three recognised marked-quartic invariants to a verified equation over Q.
//
//   magma -b DISC:=35 LEV:=1 WS:=1 UM:=7 AA:=1939/952200 BB:=34500/76729 CC0:=1939/5000 \
//          hypfinish3.m < /dev/null
//
// WHAT IS LEFT AFTER THE INVARIANTS, AND WHY IT IS NOT SEARCHED FOR.
// A, B, C fix the quartic only up to t -> lambda t and rescaling, so the curve is
//     y^2 = kappa * gnorm(mu x^2),   gnorm(t) = 1 + t + A t^2 + A^2 B t^3 + A^3 B^2 C t^4,
// with kappa and mu rational and each mattering only modulo squares.  A complex-analytic
// computation cannot see either: the coordinate it chose differs from a rational one by a
// complex scalar, which moves kappa and mu by complex squares.  Both are arithmetic and both
// are SOLVED rather than tried:
//   * kappa is read off the genus-1 Atkin-Lehner quotient, whose equation is Y^2 = kappa gnorm(t)
//     -- mu disappears there under t -> t/mu -- so each prime gives the Legendre symbol
//     (kappa|p) and one linear equation over F_2;
//   * mu is then read off the curve itself the same way, kappa being known.
// Both systems are heavily overdetermined and the solution must be unique; a kernel means the
// primes do not separate the candidates and the cap has to be raised.
AttachSpec("/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/ShimuraQuotients.spec");
SetColumns(0);
D := StringToInteger(DISC);
N := (assigned LEV) select StringToInteger(LEV) else 1;
W := { StringToInteger(x) : x in Split(WS, ",") };
um := StringToInteger(UM);
A := StringToRational(AA); B := StringToRational(BB); Cv := StringToRational(CC0);
CAP := StringToInteger(assigned PCAP select PCAP else "150");
PQ<x> := PolynomialRing(Rationals());

gnorm := 1 + x + A*x^2 + A^2*B*x^3 + A^3*B^2*Cv*x^4;
printf "=== X_0^%o(%o)/%o, u = w_%o\n", D, N, Sort(SetToSequence(W)), um;
printf "gnorm(t) = %o\n", gnorm;
error if Degree(gnorm) ne 4 or Discriminant(gnorm) eq 0, "the quartic is degenerate";

// the Atkin-Lehner overgroup, whose quotient has genus 1
WU := { (a*um) div GCD(a,um)^2 : a in W } join W;
gE := GenusShimuraCurveQuotient(D, N, WU);
gX := GenusShimuraCurveQuotient(D, N, W);
printf "X genus %o; X/w_%o has group %o and genus %o\n", gX, um, Sort(SetToSequence(WU)), gE;
error if gX ne 3 or gE ne 1, "expected genus 3 over a genus-1 Atkin-Lehner quotient";

XQ := CreateShimuraQuot(D, N, W);      XQ`g := gX;
EQ := CreateShimuraQuot(D, N, WU);     EQ`g := gE;
good := [ p : p in PrimesUpTo(CAP) | p ne 2 and (D*N) mod p ne 0 ];
tX := AssociativeArray(); tE := AssociativeArray();
for p in good do
    tX[p] := [ ComputePointsViaTrace(XQ, p, d) : d in [1,2] ];
    tE[p] := [ ComputePointsViaTrace(EQ, p, d) : d in [1,2] ];
end for;
printf "trace-formula counts taken at %o good odd primes below %o\n", #good, CAP;

// support of the two unknown numbers: bad primes of the quartic and of 2DN
function Supp(poly)
    cs := Coefficients(poly);
    dd := LCM([ Denominator(c) : c in cs ]);
    ip := Integers() ! (dd^4 * Discriminant(poly));
    S := Set(PrimeDivisors(2*D*N));
    for q in Factorisation(AbsoluteValue(ip)) do Include(~S, q[1]); end for;
    for c in cs do
        n := Numerator(c); d := Denominator(c);
        if n ne 0 then for q in Factorisation(AbsoluteValue(n)) do Include(~S, q[1]); end for; end if;
        for q in Factorisation(d) do Include(~S, q[1]); end for;
    end for;
    return Sort(SetToSequence(S));
end function;
supp := Supp(gnorm);
printf "the two numbers can only be supported on %o primes: %o\n", #supp, supp;

// Scaling y^2 = f(x) by the square of a common denominator is an isomorphism over Q (y scales
// with it), so the integral model has the same point counts; a prime dividing that denominator
// is reported as bad rather than counted on a model it does not belong to.
function CountQ(poly, p, d)
    den := LCM([ Denominator(c) : c in Coefficients(poly) ] cat [1]);
    if den mod p eq 0 then return -1; end if;
    pint := PolynomialRing(Rationals()) ! (den^2 * poly);
    F := GF(p^d);
    gg := ChangeRing(pint, F);
    if Degree(gg) lt Degree(pint) - 1 or Degree(gg) lt 3 or Discriminant(gg) eq 0 then
        return -1;
    end if;
    return #Points(HyperellipticCurve(gg));
end function;

// --- kappa, from the genus-1 quotient ---
rows := []; rhs := []; nus := 0;
for p in good do
    n1 := CountQ(gnorm, p, 1);
    if n1 lt 0 then continue; end if;
    t1 := tE[p][1];
    if n1 eq t1 and 2*(p+1) - n1 eq t1 then continue; end if;
    if n1 eq t1 then b := 0;
    elif 2*(p+1) - n1 eq t1 then b := 1;
    else nus +:= 1; continue; end if;
    Append(~rows, [ GF(2) | (LegendreSymbol(q, p) eq 1 select 0 else 1) : q in supp ]
                  cat [ GF(2) | (p mod 4 eq 1 select 0 else 1) ]);
    Append(~rhs, GF(2) ! b);
end for;
printf "kappa: %o primes give a symbol, %o unusable, %o unknowns\n", #rows, nus, #supp+1;
error if #rows lt #supp + 4, "too few primes for kappa";
Mk := Matrix(GF(2), #rows, #supp+1, &cat rows);
okk, vk := IsConsistent(Transpose(Mk), Vector(GF(2), rhs));
error if not okk, "the symbols for kappa are inconsistent: the quartic is not a twist of X/u";
printf "kappa kernel dimension %o\n", #supp + 1 - Rank(Mk);
error if Rank(Mk) ne #supp+1, "kappa is not determined uniquely; raise the prime cap";
ksgn := (vk[#supp+1] eq GF(2)!1) select -1 else 1;
kap := ksgn * (&*([Integers() | 1] cat [ supp[i] : i in [1..#supp] | vk[i] eq GF(2)!1 ]));
printf "kappa = %o\n", kap;
Ecur := HyperellipticCurve(kap*gnorm);
nE := 0;
for p in good do for d in [1,2] do
    c := CountQ(kap*gnorm, p, d);
    if c lt 0 then continue; end if;
    error if c ne tE[p][d], Sprintf("X/u count disagrees at p=%o d=%o: %o against %o", p, d, c, tE[p][d]);
    nE +:= 1;
end for; end for;
printf "X_0^%o(%o)/w_%o VERIFIED at %o counts: %o\n", D, N, um, nE, Ecur;

// --- mu, from the curve itself ---
rows := []; rhs := []; nus := 0;
for p in good do
    f0 := kap*Evaluate(gnorm, x^2);           // mu a square mod p
    nsq := 1; while LegendreSymbol(nsq, p) eq 1 do nsq +:= 1; end while;
    f1 := kap*Evaluate(gnorm, nsq*x^2);       // mu a non-square mod p
    n0 := CountQ(f0, p, 1); n1 := CountQ(f1, p, 1);
    if n0 lt 0 or n1 lt 0 or n0 eq n1 then continue; end if;
    t1 := tX[p][1];
    if n0 eq t1 then b := 0; elif n1 eq t1 then b := 1; else nus +:= 1; continue; end if;
    Append(~rows, [ GF(2) | (LegendreSymbol(q, p) eq 1 select 0 else 1) : q in supp ]
                  cat [ GF(2) | (p mod 4 eq 1 select 0 else 1) ]);
    Append(~rhs, GF(2) ! b);
end for;
printf "mu: %o primes give a symbol, %o unusable, %o unknowns\n", #rows, nus, #supp+1;
error if #rows lt #supp + 4, "too few primes for mu";
Mm := Matrix(GF(2), #rows, #supp+1, &cat rows);
okm, vm := IsConsistent(Transpose(Mm), Vector(GF(2), rhs));
error if not okm, "the symbols for mu are inconsistent";
printf "mu kernel dimension %o\n", #supp + 1 - Rank(Mm);
error if Rank(Mm) ne #supp+1, "mu is not determined uniquely; raise the prime cap";
msgn := (vm[#supp+1] eq GF(2)!1) select -1 else 1;
mu := msgn * (&*([Integers() | 1] cat [ supp[i] : i in [1..#supp] | vm[i] eq GF(2)!1 ]));
printf "mu = %o\n", mu;

foct := kap*Evaluate(gnorm, mu*x^2);
dd := LCM([ Denominator(c) : c in Coefficients(foct) ]);
foct := dd^2*foct;
Cfin := HyperellipticCurve(foct);
try Cfin := ReducedMinimalWeierstrassModel(Cfin); catch e
    try Cfin := ReducedModel(SimplifiedModel(Cfin)); catch e2 ; end try;
end try;
printf "\nMODEL X_0^%o(%o)/%o :\n   %o\n", D, N, Sort(SetToSequence(W)), Cfin;
printf "genus %o (formula says %o)\n", Genus(Cfin), gX;
error if Genus(Cfin) ne gX, "genus mismatch";

fq, hq := HyperellipticPolynomials(Cfin);
gg := fq + hq^2/4;
nchk := 0; nskip := 0;
for p in good do for d in [1,2] do
    c := CountQ(gg, p, d);
    if c lt 0 then nskip +:= 1; continue; end if;
    error if c ne tX[p][d],
       Sprintf("POINT COUNT DISAGREES at p = %o, degree %o: %o against %o", p, d, c, tX[p][d]);
    nchk +:= 1;
end for; end for;
printf "VERIFIED %o point counts over F_p and F_p^2 against the trace formula (%o skipped "
     * "for bad reduction), at %o primes below %o\n", nchk, nskip, #good, CAP;

fn := Sprintf("/Users/assaferan/Documents/GitHub/ShimuraCurveALQuotients/worktrees/campaign/data/models/models_%o_%o.m", D, N);
if FileExists(fn) then
    models := eval (Read(fn) cat "\nreturn models;");
    key := Sort(SetToSequence(W));
    if IsDefined(models, key) and #models[key] ne 0 and Type(models[key][1][2]) ne MonStgElt then
        e := models[key][1];
        Cold := HyperellipticCurve(e[2] + (#e ge 3 select e[3]^2/4 else 0));
        iso := IsIsomorphic(Cold, Cfin);
        printf "STORED MODEL: %o\n", Cold;
        printf "ISOMORPHIC OVER Q TO THE STORED MODEL: %o\n", iso;
    else
        printf "no stored hyperelliptic model at this key\n";
    end if;
else
    printf "no model file for this base\n";
end if;
printf "HYPFINISH3 DONE\n";
quit;
