// tests/GuoYangQuotients_26_3.m -- complete oracle for X_0^26(3), derived by hand.
//
// ⚠ WHY BY HAND: Guo-Yang present this base as a CRV PAIR (z^2 = -8x^2-3, y^2 = x^6-2x^4+9x^2+8),
// and CurveQuotient cannot act on a curve in a weighted projective ambient -- Magma models that as
// a toric variety and IdentityMap returns a TorMap, so AutomorphismGroup dies
// (Magma-Maths/Magma#123). Same situation as 10_19 and 10_13.
//
// ⚠ AND THERE IS NO CM TABLE FOR THIS BASE to fall back on: 26_3 is one of the two Guo-Yang CM
// tables never added as an offline test, because of the s <-> s~ swap at discs -267 and -708. So
// the fixed-point argument that settled 10_13's labelling is not available here. It is not needed:
// every label matches on the first try.
//
// Guo-Yang publish w_2 = (-x,-y,-z), w_3 = (x,-y,-z), w_26 = (x,-y,z). Composing (indices multiply
// and divide out the square of the gcd, among the divisors of 78):
//     w_6  = w_2 w_3  = (-x,  y,  z)      w_13 = w_2 w_26 = (-x,  y, -z)
//     w_78 = w_3 w_26 = ( x,  y, -z)      w_39 = w_6 w_26 = (-x, -y,  z)
// The action is DIAGONAL and both polynomials are EVEN in x, so with u = x^2 they are F(u), G(u),
// and each quotient is the field of invariant monomials. G is LINEAR, so z^2 = G(u) inverts and
// the quotients in which z survives collapse to plain hyperelliptic models; the conic
// c^2 = u*G(u) = u(-8u-3) has the rational point (0,0), so it parametrises by c = t*u for the rest.

gys_Pu<u> := PolynomialRing(Rationals());
gys_F := u^3 - 2*u^2 + 9*u + 8;
gys_G := -8*u - 3;
gys_Px<x> := PolynomialRing(Rationals());
gys_f := Evaluate(gys_F, x^2);
gys_g := Evaluate(gys_G, x^2);
gys_K<X> := FunctionField(Rationals());
error if Evaluate(gys_F, X^2) ne Evaluate(gys_f, X), "26_3: F(u) does not reproduce f";
error if Evaluate(gys_G, X^2) ne Evaluate(gys_g, X), "26_3: G(u) does not reproduce g";

gys_Pz<z> := PolynomialRing(Rationals());
gys_sub := -(z^2 + 3)/8;
gys_Kt<t> := FunctionField(Rationals());
gys_uu := -3/(t^2 + 8);
gys_Pt<T> := PolynomialRing(Rationals());
// ⚠ Clearing denominators leaves SQUARE factors, and y^2 = Q^2*R is singular as written; strip
// factors of even multiplicity (Y = y/Q), keeping the leading constant, which fixes the twist.
function gys_sqfree(e)
    error if Denominator(e) ne 1, "26_3: parametrised expression is not a polynomial";
    P := Evaluate(gys_Pt ! Numerator(e), T);
    sq := gys_Pt ! 1;
    for fe in Factorisation(P) do sq *:= fe[1]^(fe[2] div 2); end for;
    return P div (sq^2);
end function;

gys_oracle := [*
  <[1,78], HyperellipticCurve(gys_Pz ! gys_f),                         "(x, y)">,
  <[1,26], HyperellipticCurve(gys_Pz ! gys_g),                         "(x, z)">,
  <[1,3],  HyperellipticCurve(gys_Pz ! (gys_f*gys_g)),                 "(x, yz)">,
  <[1,6],  HyperellipticCurve(Evaluate(gys_F, gys_sub)),               "(u, y) via z">,
  <[1,39], HyperellipticCurve(gys_sub*Evaluate(gys_F, gys_sub)),       "(u, xy) via z">,
  <[1,2],  HyperellipticCurve(gys_sqfree((t^2+8)^4*(gys_uu*Evaluate(gys_F, gys_uu)))),
                                                                       "(u, xy, xz)">,
  <[1,13], HyperellipticCurve(gys_sqfree((t^2+8)^4*Evaluate(gys_F, gys_uu))),
                                                                       "(u, y, xz)">
*];

// ⚠ A STORED ENTRY MAY BE <genus, f, h>, MEANING y^2 + h*y = f (see GuoYangQuotientOracle.m).
function gys_model_curve(e)
    if (#e ge 3) and (Type(e[3]) eq RngUPolElt) and (e[3] ne 0) then
        return HyperellipticCurve(e[2], e[3]);
    end if;
    return HyperellipticCurve(e[2]);
end function;

gys_models := eval (Read("data/models/models_26_3.m") cat "\nreturn models;");
gys_n := 0; gys_empty := 0;
for gys_o in gys_oracle do
    gys_lab, gys_Cq, gys_gen := Explode(gys_o);
    gys_ok, gys_es := IsDefined(gys_models, [Integers()| s : s in gys_lab]);
    if (not gys_ok) or (#gys_es eq 0) then gys_empty +:= 1; continue; end if;
    for gys_e in gys_es do
        if Type(gys_e[2]) eq MonStgElt then continue; end if;
        gys_Cs := gys_model_curve(gys_e);
        error if Genus(gys_Cs) ne Genus(gys_Cq),
            Sprintf("X0^26(3) W=%o: our genus %o vs Guo-Yang's %o -- wrong object",
                    gys_lab, Genus(gys_Cs), Genus(gys_Cq));
        gys_r := (Genus(gys_Cq) eq 0)
                 select HasRationalPoint(Conic(gys_Cs)) eq HasRationalPoint(Conic(gys_Cq))
                 else IsIsomorphic(gys_Cs, gys_Cq);
        error if not gys_r,
            Sprintf("X0^26(3) W=%o: our stored curve is NOT isomorphic to Guo-Yang's quotient by "
                    * "the invariants %o", gys_lab, gys_gen);
        gys_n +:= 1;
    end for;
end for;

// ⚠ [1,2] and [1,13] are the two EMPTY keys at this base. When EquationsByRebase fills them, this
// oracle already holds the curves they must equal -- so they get checked the moment they appear
// rather than being accepted because they look plausible.
error if gys_n lt 7,
    Sprintf("X0^26(3): expected at least 7 quotient comparisons, made %o (%o empty)",
            gys_n, gys_empty);
printf " ok (X0^26(3): %o quotient(s) checked against Guo-Yang's curve + involutions, "
       * "%o key(s) still empty)\n", gys_n, gys_empty;
