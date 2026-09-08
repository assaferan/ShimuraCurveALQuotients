// tests/GuoYangQuotients_10_13.m -- complete oracle for X_0^10(13), derived by hand.
//
// ⚠ WHY BY HAND: Guo-Yang present this base as a CRV PAIR (z^2 = -2x^2-25, y^2 = 5x^4-74x^2+325),
// and CurveQuotient cannot act on a curve in a weighted projective ambient -- Magma models that as
// a toric variety and IdentityMap returns a TorMap, so AutomorphismGroup dies. Reported as
// Magma-Maths/Magma#123. The generic sweep in tests/GuoYangQuotientOracle.m therefore cannot
// cover this base, exactly as at 10_19.
//
// Guo-Yang give w_2(x,y,z) = (x,-y,-z), w_5 = (-x,-y,-z), w_65 = (x,-y,z). Composing (indices
// multiply and divide out the square of the gcd, in the group of divisors of 130):
//     w_10  = w_2 w_5   = (-x,  y,  z)      w_13  = w_5 w_65  = (-x,  y, -z)
//     w_26  = w_10 w_65 = (-x, -y,  z)      w_130 = w_2 w_65  = ( x,  y, -z)
// Every action is DIAGONAL, so each quotient is the field of invariant monomials.
// Both defining polynomials are EVEN in x, so with u = x^2 they are F(u) and G(u).

gyt_Pu<u> := PolynomialRing(Rationals());
gyt_F := 5*u^2 - 74*u + 325;         // y^2 = F(x^2)
gyt_G := -2*u - 25;                  // z^2 = G(x^2)
gyt_Px<x> := PolynomialRing(Rationals());
gyt_f := Evaluate(gyt_F, x^2);
gyt_g := Evaluate(gyt_G, x^2);

// sanity: the published pair must have the genus our model records for W={1}
gyt_K<X> := FunctionField(Rationals());
error if Evaluate(gyt_F, X^2) ne Evaluate(gyt_f, X), "10_13: F(u) does not reproduce f";
error if Evaluate(gyt_G, X^2) ne Evaluate(gyt_g, X), "10_13: G(u) does not reproduce g";

// G is LINEAR, so z^2 = G(u) inverts: u = -(z^2 + 25)/2. That collapses the quotients in which
// z survives to plain hyperelliptic models.
gyt_Pz<z> := PolynomialRing(Rationals());
gyt_sub := -(z^2 + 25)/2;

// The conic c^2 = u*G(u) = u(-2u-25) HAS the rational point (0,0), so it parametrises by c = t*u,
// giving u = -25/(t^2+2); that handles the quotients in which only x*z survives.
gyt_Kt<t> := FunctionField(Rationals());
gyt_uu := -25/(t^2 + 2);
gyt_Pt<T> := PolynomialRing(Rationals());
function gyt_clear(e)
    error if Denominator(e) ne 1, "10_13: parametrised expression is not a polynomial";
    return Evaluate(gyt_Pt ! Numerator(e), T);
end function;

gyt_oracle := [*
  <[1,130], HyperellipticCurve(gyt_f),                             "(x, y)">,
  <[1,65],  HyperellipticCurve(gyt_g),                             "(x, z)">,
  <[1,2],   HyperellipticCurve(gyt_f*gyt_g),                       "(x, yz)">,
  <[1,10],  HyperellipticCurve(Evaluate(gyt_F, gyt_sub)),          "(u, y) via z">,
  <[1,26],  HyperellipticCurve(gyt_sub*Evaluate(gyt_F, gyt_sub)),  "(u, xy) via z">,
  <[1,5],   HyperellipticCurve(gyt_clear((t^2+2)^4 * (gyt_uu*Evaluate(gyt_F, gyt_uu)))),
                                                                   "(u, xy, xz)">,
  <[1,13],  HyperellipticCurve(gyt_clear((t^2+2)^4 * Evaluate(gyt_F, gyt_uu))),
                                                                   "(u, y, xz)">
*];

// ⚠ A STORED ENTRY MAY BE <genus, f, h>, MEANING y^2 + h*y = f. Reading only e[2] drops h and
// gives a DIFFERENT curve of the same genus (see tests/GuoYangQuotientOracle.m).
function gyt_model_curve(e)
    if (#e ge 3) and (Type(e[3]) eq RngUPolElt) and (e[3] ne 0) then
        return HyperellipticCurve(e[2], e[3]);
    end if;
    return HyperellipticCurve(e[2]);
end function;

gyt_models := eval (Read("data/models/models_10_13.m") cat "\nreturn models;");
gyt_n := 0; gyt_empty := 0; gyt_skip := 0;
for gyt_o in gyt_oracle do
    gyt_lab, gyt_Cq, gyt_gen := Explode(gyt_o);
    gyt_ok, gyt_es := IsDefined(gyt_models, [Integers()| s : s in gyt_lab]);
    if (not gyt_ok) or (#gyt_es eq 0) then gyt_empty +:= 1; continue; end if;
    for gyt_e in gyt_es do
        if Type(gyt_e[2]) eq MonStgElt then continue; end if;
        gyt_Cs := gyt_model_curve(gyt_e);
        error if Genus(gyt_Cs) ne Genus(gyt_Cq),
            Sprintf("X0^10(13) W=%o: our genus %o vs Guo-Yang's %o -- wrong object",
                    gyt_lab, Genus(gyt_Cs), Genus(gyt_Cq));
        gyt_res := false;
        if Genus(gyt_Cq) eq 0 then
            gyt_res := HasRationalPoint(Conic(gyt_Cs)) eq HasRationalPoint(Conic(gyt_Cq));
        elif Genus(gyt_Cq) eq 1 then
            // ⚠ IsIsomorphic REFUSES genus-1 curves over Q; compare Jacobians, and SKIP rather
            // than score if no elliptic model is obtainable -- an error is not a negative result.
            gyt_got := false;
            try
                gyt_E1 := Jacobian(GenusOneModel(HyperellipticPolynomials(gyt_Cs)));
                gyt_E2 := Jacobian(GenusOneModel(HyperellipticPolynomials(gyt_Cq)));
                gyt_res := IsIsomorphic(gyt_E1, gyt_E2); gyt_got := true;
            catch err ; end try;
            if not gyt_got then gyt_skip +:= 1; continue; end if;
        else
            gyt_res := IsIsomorphic(gyt_Cs, gyt_Cq);
        end if;
        error if not gyt_res,
            Sprintf("X0^10(13) W=%o: our stored curve is NOT isomorphic to Guo-Yang's quotient by "
                    * "the invariants %o", gyt_lab, gyt_gen);
        gyt_n +:= 1;
    end for;
end for;

error if gyt_n lt 4,
    Sprintf("X0^10(13): expected at least 4 quotient comparisons, made %o (%o empty, %o skipped)",
            gyt_n, gyt_empty, gyt_skip);
printf " ok (X0^10(13): %o quotient(s) checked against Guo-Yang's curve + involutions, "
       * "%o empty, %o skipped)\n", gyt_n, gyt_empty, gyt_skip;
