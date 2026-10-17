// Jacobian of the genus-1 quotient y^2 = quartic, via the classical invariants of the quartic:
//   I = 12ae - 3bd + c^2,  J = 72ace + 9bcd - 27ad^2 - 27eb^2 - 2c^3,  Jac: Y^2 = X^3 - 27 I X - 27 J.
P<x> := PolynomialRing(Rationals());
q := -25447851/236421376*x^4 - 573417/476656*x^3 + 18083/7688*x^2 + 13*x + 13;
a := Coefficient(q,4); b := Coefficient(q,3); c := Coefficient(q,2); d := Coefficient(q,1); e := Coefficient(q,0);
I4 := 12*a*e - 3*b*d + c^2;
J4 := 72*a*c*e + 9*b*c*d - 27*a*d^2 - 27*e*b^2 - 2*c^3;
E := MinimalModel(EllipticCurve([0, 0, 0, -27*I4, -27*J4]));
printf "Jacobian of X_0^39(11)/w_3 : %o\n", E;
printf "conductor %o    (D*N = 429 = 3*11*13)\n", Conductor(E);
printf "j-invariant %o\n", jInvariant(E);
try printf "Cremona label: %o\n", CremonaReference(E); catch err printf "not in the local database\n"; end try;
printf "\nisogeny classes of conductor 429, matched by Frobenius traces:\n";
seen := {};
for F in EllipticCurves(CremonaDatabase(), 429) do
    r := CremonaReference(F); lbl := r[1..4];
    if lbl in seen then continue; end if; Include(~seen, lbl);
    same := &and[ TraceOfFrobenius(F,p) eq TraceOfFrobenius(E,p) : p in [2,5,7,17,19,23] ];
    printf "   %-8o traces match ours: %-6o isomorphic: %o\n", lbl, same, IsIsomorphic(F,E);
end for;
print "CERT4 DONE";
