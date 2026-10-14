// An Atkin-Lehner quotient X_0(D,N)/W is the fibre product over the star line of its index-2
// Atkin-Lehner double covers y^2 = f(t).  This script rebuilds such a curve from the stored
// equations of those covers, computes its canonical model (a plane quartic when the curve has
// genus 3 and is not hyperelliptic) and decides hyperellipticity two ways: Magma's IsHyperelliptic
// on the curve, and the dimension of the quadratic relations among the canonical differentials (a
// hyperelliptic genus-3 curve has canonical image a conic, so a quadratic relation; a
// non-hyperelliptic one has none and its canonical image is the quartic itself).  Point counts
// over F_p and F_{p^2} are compared with the Eichler-Selberg trace formula, which knows nothing
// about the construction; a disagreement at a good prime stops the script.
//
//     magma -b D:=6 N:=23 W:=23 tools/planequartic.m          (X_0^6(23)/w_23)
//     magma -b D:=34 N:=3 W:=2 tools/planequartic.m           (X_0^34(3)/w_2)
//     magma -b D:=46 N:=3 W:=3 tools/planequartic.m           (X_0^46(3)/w_3)
//
// The three curves named first are the ones recorded as non-hyperelliptic in data/models/PROVENANCE.md;
// this is the computation behind that record, and tests/PlaneQuartics.m
// is the test that keeps it true.
AttachSpec("ShimuraQuotients.spec");
SetColumns(0);
D := StringToInteger(D); N := StringToInteger(N); w := StringToInteger(W);
DN := D*N;
models := eval (Read(Sprintf("data/models/models_%o_%o.m", D, N)) cat "\nreturn models;");
Wset := {Integers()| 1, w};
full := {Integers()| d : d in Divisors(DN) | GCD(d, DN div d) eq 1};
error if Wset notsubset full, Sprintf("w_%o is not an Atkin-Lehner involution of X_0(%o,%o)", w, D, N);
// the double covers of the star curve above W: the Atkin-Lehner groups of index 2 containing W
sups := AtkinLehnerDoubleCoversOver(Wset, D, N);
P<x> := PolynomialRing(Rationals());
avail := {}; poly := AssociativeArray();
for U in sups do
    k := Sort(SetToSequence(U));
    if not IsDefined(models, k) or #models[k] eq 0 or Type(models[k][1][2]) eq MonStgElt then continue; end if;
    e := models[k][1];
    // y^2 + h y = f is (y + h/2)^2 = f + h^2/4
    poly[U] := e[2] + (#e ge 3 select e[3]^2/4 else 0);
    Include(~avail, U);
    printf "  above W: %o, genus %o\n", k, e[1];
end for;
// as many of them as the index needs, meeting in exactly W
ok, gens := FibreProductGenerators(Wset, avail, D, N);
error if not ok,
    Sprintf("the %o stored double covers above w_%o do not generate X_0^%o(%o)/w_%o (index %o)",
            #avail, w, D, N, w, #full div #Wset);
fs := [poly[U] : U in gens];
printf "X_0^%o(%o)/w_%o: fibre product of %o\n", D, N, w, [Sort(SetToSequence(U)) : U in gens];
A := AffineSpace(Rationals(), 1 + #fs);
t := A.1;
C := Curve(A, [A.(1+i)^2 - Evaluate(fs[i], t) : i in [1..#fs]]);
Cp := ProjectiveClosure(C);
g := Genus(Cp);
printf "genus %o\n", g;
// point counts against the trace formula, at the first good primes: a prime where the tower
// cannot be built or the genus drops is bad reduction of this model and is replaced by the next
curves := GetHyperellipticCandidates();
X := rep{Y : Y in curves | Y`D eq D and Y`N eq N and Y`W eq Wset};
error if X`g ne g, Sprintf("the genus formula gives %o, the fibre product has genus %o", X`g, g);
checked := 0; p := 2; tried := 0;
while checked lt 4 and tried lt 40 do
    p := NextPrime(p); tried +:= 1;
    if DN mod p eq 0 then continue; end if;
    Kp := RationalFunctionField(GF(p)); L := Kp; built := true;
    try
        for f in fs do
            R<Y> := PolynomialRing(L);
            L := FunctionField(Y^2 - L!Evaluate(PolynomialRing(GF(p))!f, Kp.1));
        end for;
    catch err built := false;
    end try;
    if not built or Genus(L) ne g then
        printf "  p = %o: bad reduction of this model, skipped\n", p;
        continue;
    end if;
    cnt := [&+[e * #Places(L, e) : e in Divisors(d)] : d in [1..2]];
    exp := [ComputePointsViaTrace(X, p, d) : d in [1..2]];
    error if cnt ne exp,
        Sprintf("p = %o: points over F_p, F_p^2 = %o, the trace formula gives %o", p, cnt, exp);
    printf "  p = %o: points over F_p, F_p^2 = %o, trace formula %o ok\n", p, cnt, exp;
    checked +:= 1;
end while;
error if checked lt 2, Sprintf("only %o good primes among the first 40", checked);
// the canonical model and hyperellipticity
hyp := IsHyperelliptic(Cp);
printf "IsHyperelliptic: %o\n", hyp;
phi := CanonicalMap(Cp);
Img := CanonicalImage(Cp, phi);
printf "canonical image in P^%o: degree %o\n", Dimension(Ambient(Img)), Degree(Img);
if g eq 3 then
    printf "  %o\n", Dimension(Ambient(Img)) eq 2 and Degree(Img) eq 4 select "a plane quartic" else "a conic (hyperelliptic)";
    if Dimension(Ambient(Img)) eq 2 and Degree(Img) eq 4 then
        printf "  nonsingular: %o; equation: %o\n", IsNonsingular(Img), DefiningPolynomial(Img);
    end if;
end if;
