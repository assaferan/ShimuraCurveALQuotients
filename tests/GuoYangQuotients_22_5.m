// tests/GuoYangQuotients_22_5.m -- a COMPLETE oracle for X_0^22(5), the companion to
// tests/GuoYangQuotients_10_19.m.
//
// Guo-Yang print only the FULL curve for 22_5 (no worked example, no quotient equations), so
// none of our fifteen cover keys looked checkable against them. But they also print the
// INVOLUTIONS, and every quotient follows from those.
//
// ⚠ UNLIKE 10_19, THE ACTION IS NOT DIAGONAL: w_2 is a Mobius map. So the quotients come from
// INVARIANT FUNCTIONS, not sign patterns. Guo-Yang give
//     w_2 :(x,y) -> ( 1/x,  y/x^6),  w_5 :(x,y) -> (-1/x, -y/x^6),  w_110:(x,y) -> (x,-y);
// composing, w_10 = w_2 w_5 : (x,y) -> (-x,-y) and w_11 = w_10 w_110 : (x,y) -> (-x, y), and then
//     w_22 = (-1/x,  y/x^6),  w_55 = ( 1/x, -y/x^6).
// Each is CHECKED below to preserve y^2 = f(x) before being used.
//
// f is EVEN and PALINDROMIC, which is what makes this clean: with T = x^2 + x^-2 we get
// f/x^6 = Q(T), and with u = x^2 we get f = F(u). Both identities are VERIFIED here rather than
// asserted, since the whole oracle rests on them.

gyv_Px<x> := PolynomialRing(Rationals());
gyv_f := -11*x^12 - 80*x^10 - 240*x^8 - 362*x^6 - 240*x^4 - 80*x^2 - 11;
error if Genus(HyperellipticCurve(gyv_f)) ne 5, "X0^22(5): Guo-Yang's curve should have genus 5";

gyv_K<X> := FunctionField(Rationals());
// <label, image of x, factor on y, sign on y>
gyv_acts := [* <2, 1/X, 1/X^6, 1>, <5, -1/X, -1/X^6, 1>, <110, X, 1, -1>,
               <10, -X, 1, -1>, <11, -X, 1, 1>, <22, -1/X, 1/X^6, 1>, <55, 1/X, -1/X^6, 1> *];
for gyv_a in gyv_acts do
    gyv_m, gyv_xi, gyv_yf, gyv_s := Explode(gyv_a);
    error if (gyv_s*gyv_yf)^2 * Evaluate(gyv_f, X) ne Evaluate(gyv_f, gyv_xi),
        Sprintf("X0^22(5): w_%o does not preserve Guo-Yang's curve -- the group law is wrong", gyv_m);
end for;

gyv_PT<T> := PolynomialRing(Rationals());
gyv_Q := -11*(T^3 - 3*T) - 80*(T^2 - 2) - 240*T - 362;
error if Evaluate(gyv_Q, X^2 + 1/X^2) ne Evaluate(gyv_f, X)/X^6,
    "X0^22(5): Q(T) does not reproduce f/x^6";
gyv_PU<u> := PolynomialRing(Rationals());
gyv_F := -11*u^6 - 80*u^5 - 240*u^4 - 362*u^3 - 240*u^2 - 80*u - 11;
error if Evaluate(gyv_F, X^2) ne Evaluate(gyv_f, X), "X0^22(5): F(u) does not reproduce f";

gyv_Pv<v> := PolynomialRing(Rationals());
gyv_oracle := [*
  <[1,11], HyperellipticCurve(Evaluate(gyv_F, v)),                    "(u = x^2, y)">,
  <[1,10], HyperellipticCurve(v*Evaluate(gyv_F, v)),                  "(u, xy)">,
  <[1,2],  HyperellipticCurve(Evaluate(gyv_Q, v^2-2)),                "(x+1/x, y/x^3)">,
  <[1,5],  HyperellipticCurve(Evaluate(gyv_Q, v^2+2)),                "(x-1/x, y/x^3)">,
  <[1,22], HyperellipticCurve(Evaluate(gyv_Q, v^2+2)*(v^2+4)),        "(x-1/x, (y/x^3)(x+1/x))">,
  <[1,55], HyperellipticCurve(Evaluate(gyv_Q, v^2-2)*(v^2-4)),        "(x+1/x, (y/x^3)(x-1/x))">
*];

gyv_models := eval (Read("data/models/models_22_5.m") cat "\nreturn models;");
gyv_n := 0; gyv_empty := 0;
for gyv_o in gyv_oracle do
    gyv_lab, gyv_Cq, gyv_gen := Explode(gyv_o);
    gyv_ok, gyv_es := IsDefined(gyv_models, [Integers()| t : t in gyv_lab]);
    if (not gyv_ok) or (#gyv_es eq 0) then gyv_empty +:= 1; continue; end if;
    for gyv_e in gyv_es do
        if Type(gyv_e[2]) eq MonStgElt then continue; end if;
        gyv_Cs := HyperellipticCurve(gyv_e[2]);
        error if Genus(gyv_Cs) ne Genus(gyv_Cq),
            Sprintf("X0^22(5) W=%o: our genus %o vs Guo-Yang's %o -- wrong object",
                    gyv_lab, Genus(gyv_Cs), Genus(gyv_Cq));
        error if not IsIsomorphic(gyv_Cs, gyv_Cq),
            Sprintf("X0^22(5) W=%o: our stored curve is NOT isomorphic to Guo-Yang's quotient by "
                    * "the invariants %o", gyv_lab, gyv_gen);
        gyv_n +:= 1;
    end for;
end for;

// [1,110] is the bare x-line -- a P^1 over Q, which is exactly WHY 22_5 is hyperelliptic over Q
// (Guo-Yang, Remark 38). Checked as "genus 0 with a rational point", not by isomorphism.
gyv_ok110, gyv_es110 := IsDefined(gyv_models, [Integers()|1,110]);
error if not gyv_ok110, "X0^22(5): no [1,110] key";
for gyv_e in gyv_es110 do
    if Type(gyv_e[2]) eq MonStgElt then continue; end if;
    gyv_C := HyperellipticCurve(gyv_e[2]);
    error if Genus(gyv_C) ne 0, "X0^22(5) W=[1,110]: expected genus 0";
    error if not HasRationalPoint(Conic(gyv_C)),
        "X0^22(5) W=[1,110]: X/w_110 must be a P^1 over Q -- that is why the curve is "
        * "hyperelliptic over Q (Guo-Yang, Remark 38)";
    gyv_n +:= 1;
end for;

error if gyv_n lt 6,
    Sprintf("X0^22(5): expected at least 6 quotient comparisons, made %o", gyv_n);
printf " ok (X0^22(5): %o quotient(s) checked against Guo-Yang's curve + involutions, "
       * "%o key(s) still empty)\n", gyv_n, gyv_empty;
