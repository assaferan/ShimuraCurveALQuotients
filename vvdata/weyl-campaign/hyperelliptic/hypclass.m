// For each hyperelliptic CRV curve X_0(D,N)/W: is the hyperelliptic involution an Atkin-Lehner one
// (some index-2 AL overgroup U of W has a genus-0 quotient)?  If not, the involution is new.  For
// genus 3 the canonical image is a plane conic: does it have a rational point?
AttachSpec("ShimuraQuotients.spec");
SetColumns(0);
P<x> := PolynomialRing(Rationals());
curves := GetHyperellipticCandidates();
lines := Split(Read(LIST), "\n");
for l in lines do
    parts := Split(l, " ");
    b := parts[1]; D := StringToInteger(Split(b, "_")[1]); N := StringToInteger(Split(b, "_")[2]);
    W := {Integers()| StringToInteger(s) : s in Split(parts[2], ",")};
    X := rep{Y : Y in curves | Y`D eq D and Y`N eq N and Y`W eq W};
    DN := D*N;
    full := {Integers()| d : d in Divisors(DN) | GCD(d, DN div d) eq 1};
    ups := AtkinLehnerDoubleCoversOver(W, D, N);
    g0 := [U : U in ups | exists{Y : Y in curves | Y`D eq D and Y`N eq N and Y`W eq U and Y`g eq 0}];
    models := eval (Read(Sprintf("data/models/models_%o.m", b)) cat "\nreturn models;");
    e := models[Sort(SetToSequence(W))][1];
    shape := exists{st : st in e[3] | Regexp("y[0-9]", st)} select "fibre product" else "pair";
    status := IsEmpty(g0) select "involution NOT Atkin-Lehner" else Sprintf("AL involution, genus-0 quotient by %o", [Sort(SetToSequence(U)) : U in g0]);
    extra := "";
    if X`g eq 3 then
        fs := CRVDoubleCovers(e[3]);
        A := AffineSpace(Rationals(), #fs + 1);
        C := ProjectiveClosure(Curve(A, [A.(i+1)^2 - Evaluate(fs[i], A.1) : i in [1..#fs]]));
        try
            phi := CanonicalMap(C); R := CanonicalImage(C, phi);
            Q := Conic(R);
            extra := Sprintf("; canonical conic has rational point: %o", HasRationalPoint(Q));
        catch err extra := "; canonical image failed"; end try;
    end if;
    printf "%o %o genus %o (%o): %o%o\n", b, Sort(SetToSequence(W)), X`g, shape, status, extra;
end for;
quit;
