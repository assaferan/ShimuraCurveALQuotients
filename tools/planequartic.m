// A genus-3 Atkin-Lehner quotient whose three index-2 quotients all have genus 1 is the fibre
// product over the star line of any two of them.  This script rebuilds such a curve from the two
// stored genus-1 equations, computes its canonical model (a plane quartic when the curve is not
// hyperelliptic) and decides hyperellipticity two ways: Magma's IsHyperelliptic on the curve, and
// the dimension of the quadratic relations among the canonical differentials (a hyperelliptic
// genus-3 curve has canonical image a conic, so a quadratic relation; a non-hyperelliptic one has
// none and its canonical image is the quartic itself).  Point counts over F_p and F_{p^2} are
// compared with the Eichler-Selberg trace formula, which knows nothing about the construction.
//
//     magma -b D:=6 N:=23 W:=23 tools/planequartic.m          (X_0^6(23)/w_23)
//     magma -b D:=34 N:=3 W:=2 tools/planequartic.m           (X_0^34(3)/w_2)
//     magma -b D:=46 N:=3 W:=3 tools/planequartic.m           (X_0^46(3)/w_3)
//
// The three curves named first are the ones recorded as non-hyperelliptic in data/models/PROVENANCE.md
// and HANDOFF.md (2026-10-03); this is the computation behind that record.
AttachSpec("ShimuraQuotients.spec");
SetColumns(0);
D := StringToInteger(D); N := StringToInteger(N); w := StringToInteger(W);
DN := D*N;
models := eval (Read(Sprintf("data/models/models_%o_%o.m", D, N)) cat "\nreturn models;");
Wset := {1, w};
full := {Integers()| d : d in Divisors(DN) | GCD(d, DN div d) eq 1};
// the index-2 Atkin-Lehner groups above <w>: <w, u> for u in full
sups := {Sort(SetToSequence({x*y div GCD(x,y)^2 : x in Wset, y in {1, u}})) : u in full | u notin Wset};
P<x> := PolynomialRing(Rationals());
polys := [];
for k in sups do
    if not IsDefined(models, k) or #models[k] eq 0 or Type(models[k][1][2]) eq MonStgElt then continue; end if;
    e := models[k][1];
    g := e[1]; f := e[2]; h := #e ge 3 select e[3] else P!0;
    if g ne 1 then continue; end if;
    Append(~polys, <k, f + h^2/4>);
end for;
printf "X_0^%o(%o)/w_%o: %o genus-1 quotients with stored equations: %o\n", D, N, w, #polys, [t[1] : t in polys];
error if #polys lt 2, "need two stored genus-1 quotients";
f1 := polys[1][2]; f2 := polys[2][2];
A<t, y1, y2> := AffineSpace(Rationals(), 3);
C := Curve(A, [y1^2 - Evaluate(f1, t), y2^2 - Evaluate(f2, t)]);
Cp := ProjectiveClosure(C);
g := Genus(Cp);
printf "fibre product of %o and %o: genus %o\n", polys[1][1], polys[2][1], g;
// point counts against the trace formula
curves := GetHyperellipticCandidates();
X := rep{Y : Y in curves | Y`D eq D and Y`N eq N and Y`W eq Wset};
assert X`g eq g;
for p in [q : q in [5, 7, 11, 13] | DN mod q ne 0] do
    Kp := RationalFunctionField(GF(p)); L := Kp;
    for f in [f1, f2] do
        R<Y> := PolynomialRing(L);
        L := FunctionField(Y^2 - L!Evaluate(PolynomialRing(GF(p))!f, Kp.1));
    end for;
    cnt := [&+[e * #Places(L, e) : e in Divisors(d)] : d in [1..2]];
    exp := [ComputePointsViaTrace(X, p, d) : d in [1..2]];
    printf "  p = %o: points over F_p, F_p^2 = %o, trace formula %o %o\n", p, cnt, exp, cnt eq exp select "ok" else "MISMATCH";
end for;
// the canonical model and hyperellipticity
hyp := IsHyperelliptic(Cp);
printf "IsHyperelliptic: %o\n", hyp;
phi := CanonicalMap(Cp);
Img := CanonicalImage(Cp, phi);
printf "canonical image in P^%o: degree %o, %o\n", Dimension(Ambient(Img)), Degree(Img), Dimension(Ambient(Img)) eq 2 and Degree(Img) eq 4 select "a plane quartic" else "a conic (hyperelliptic)";
if Dimension(Ambient(Img)) eq 2 and Degree(Img) eq 4 then
    printf "  nonsingular: %o; equation: %o\n", IsNonsingular(Img), DefiningPolynomial(Img);
end if;
