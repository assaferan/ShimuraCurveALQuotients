// Fit the constant: compare the canonical heights of the known branch points with
// L(E,1) * L'(E^d,1) * sqrt|d| / Area(E), where Area is the area of the period lattice of E.
AttachSpec("ShimuraQuotients.spec"); SetColumns(0); SetDefaultRealFieldPrecision(30);
P<x> := PolynomialRing(Rationals());
for D in [35, 51] do
    models := eval (Read("data/models/models_" cat IntegerToString(D) cat "_1.m") cat " return models;");
    Ftop := models[[Integers()|1]][1][2];
    q := P![Coefficient(Ftop, 2*i) : i in [0..4]];
    den := LCM([Denominator(c) : c in Coefficients(q)]);
    C := HyperellipticCurve(q*den^2);
    E0, mp := EllipticCurve(C, Points(C : Bound := 5000)[1]);
    E, mn := MinimalModel(E0);
    u := Pullback(Inverse(mp), FunctionField(C).1);
    Z, Pl := SignDecomposition(Divisor(u));
    per := Periods(E); A := Abs(Imaginary(ComplexConjugate(per[1])*per[2]));
    printf "\n=== D=%o: E %o (Cremona %o), area of period lattice %o, real period %o ===\n", D, aInvariants(E), CremonaReference(E), A, RealPeriod(E);
    for pl in Support(Z + Pl) do
        R := RepresentativePoint(pl); K := Ring(Parent(R));
        d := Discriminant(Integers(AbsoluteField(K)));
        EK := BaseChange(E, K); PK := EK!Eltseq(mn(E0(K)!R));
        h := Height(PK);
        tw := QuadraticTwist(E, d);
        L1 := Evaluate(LSeries(E), 1); Ld1 := Evaluate(LSeries(tw), 1 : Derivative := 1);
        base := L1 * Ld1 * Sqrt(Abs(d)) / A;
        printf "d = %o: height %o, L(E,1) L'(E^d,1) sqrt|d| / Area = %o, ratio height/base = %o\n", d, h, base, h/base;
        printf "        ratio * 4pi^2 = %o, ratio / pi^2 = %o, ratio * pi^2 = %o\n", h/base*4*Pi(RealField())^2, h/base/Pi(RealField())^2, h/base*Pi(RealField())^2;
    end for;
end for;
