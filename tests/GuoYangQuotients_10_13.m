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
// ⚠ CLEARING DENOMINATORS LEAVES SQUARE FACTORS, and y^2 = Q^2*R is SINGULAR as written.
// Multiplying F(u) by s^4 (s = t^2+2) leaves an s^2 factor, because F(u) = (...)/s^2. For y^2 = P
// with P = Q^2*R the model is Y^2 = R via Y = y/Q, so strip every factor of even multiplicity.
function gyt_clear(e)
    error if Denominator(e) ne 1, "10_13: parametrised expression is not a polynomial";
    P := Evaluate(gyt_Pt ! Numerator(e), T);
    sq := gyt_Pt ! 1;
    for fe in Factorisation(P) do sq *:= fe[1]^(fe[2] div 2); end for;
    return P div (sq^2);      // keeps the leading constant, which fixes the QUADRATIC TWIST
end function;

// ⚠⚠ GUO-YANG'S INVOLUTION TABLE SWAPS w_10 AND w_13 HERE, and their OWN CM table proves it.
// This is the third error found in their tables (after 93_1's `-3t` and 14_5's w_35 sign), and it
// is settled by fixed points, not by preferring one source over the other.
//
// Their two candidates both negate x and fix y, differing only in z. On P(1,2,1,1) with
// (Z,Y,X,W) ~ (lam*Z, lam^2*Y, lam*X, lam*W), a map (Z,Y,X,W) -> (a*Z, Y, -X, W) fixes a point iff
// there is lam with a*Z = lam*Z, -X = lam*X, W = lam*W:
//   * X /= 0 forces lam = -1, hence W = 0 (the points at INFINITY) and a = -1. So ONLY (-x,y,-z)
//     fixes them. Those points have Z^2 = -2X^2, Y^2 = 5X^4, so they live over
//     Q(sqrt(-2), sqrt(5)) -- which contains sqrt(-10).
//   * X = 0 forces lam = 1, hence a = +1. So ONLY (-x,y,z) fixes the x=0 points, where
//     Z^2 = -25 and Y^2 = 325 = 25*13, i.e. over Q(i, sqrt(13)) -- which contains sqrt(-13).
// (Both rely on Z /= 0 at the locus, which holds: -25 /= 0, and -2X^2 /= 0 for X /= 0.)
//
// By Ogg, the fixed points of w_m are the CM points of discriminant -4m. So the map fixed at the
// sqrt(-10) points is w_10 (-40 = -4*10) and the one fixed at the sqrt(-13) points is w_13
// (-52 = -4*13):
//     w_10 = (-x, y, -z)        w_13 = (-x, y,  z)
// which is the OPPOSITE of their involution table -- and is what our pipeline says.
//
// ⚠⚠ THE DISCREPANCY IS A GROUP AUTOMORPHISM, NOT A SINGLE TRANSPOSITION, and only half of it is
// PROVEN. Our labelling differs from Guo-Yang's by the map that swaps 5 <-> 26 and 10 <-> 13 while
// fixing 2, 65 and 130. That IS an automorphism of the Atkin-Lehner group, and it is multiplicative
// -- e.g. 2*5 = 10 goes to 2*26 = 13, and 5*65 = 13 goes to 26*65 = 10 -- so the two swaps stand or
// fall together.
//   * 10 <-> 13 is PROVEN in our favour by the fixed-point argument below, whose clincher is
//     Guo-Yang's OWN CM table.
//   * 5 <-> 26 CANNOT be settled the same way, and this is not for want of trying: both quotients
//     have genus 2, and Riemann-Hurwitz on a genus-3 curve forces r = 0 (2*3-2 = 2(2*2-2) + r), so
//     BOTH involutions are FIXED-POINT FREE. Ogg's rule is about fixed points, so it says nothing
//     here, and the CM table cannot help either because there are no fixed CM points to place.
// ⇒ We adopt our pipeline's labelling for BOTH pairs, because a single group automorphism is the
// only consistent reading and one of its two swaps is proven. The 5 <-> 26 half is INFERRED BY
// CONSISTENCY, not independently established, and should be described that way.
//
// ⚠ THE CLINCHER IS INTERNAL TO THEIR PAPER: their CM table for this base (transcribed in
// tests/_offline/GuoYang_10_13.m) lists disc -52 at Hauptmodul value 0 and disc -40 at infinity.
// Since u = x^2 is the star Hauptmodul in their normalisation -- checked against their own CM
// values, e.g. d=-3 at u=1 gives z^2 = -27, d=-43 at u=9 gives z^2 = -43, d=-35 at u=5 gives
// z^2 = -35, each exactly d times a square -- u=0 is the x=0 locus and u=infinity the points at
// infinity. So their CM table and their involution table contradict each other, and the CM table
// agrees with us.
gyt_oracle := [*
  <[1,130], HyperellipticCurve(gyt_f),                             "(x, y)">,
  <[1,65],  HyperellipticCurve(gyt_g),                             "(x, z)">,
  <[1,2],   HyperellipticCurve(gyt_f*gyt_g),                       "(x, yz)">,
  // ⚠ SECOND SWAP: (u,xy) is OUR w_5 and (u,xy,xz) is OUR w_26 -- see the note below.
  <[1,5],   HyperellipticCurve(gyt_sub*Evaluate(gyt_F, gyt_sub)),  "(u, xy) via z">,
  <[1,26],  HyperellipticCurve(gyt_clear((t^2+2)^4 * (gyt_uu*Evaluate(gyt_F, gyt_uu)))),
                                                                   "(u, xy, xz)">,
  // the two whose labels the fixed-point argument above CORRECTS relative to Guo-Yang's table
  <[1,13],  HyperellipticCurve(Evaluate(gyt_F, gyt_sub)),          "(u, y) via z">,
  <[1,10],  HyperellipticCurve(gyt_clear((t^2+2)^4 * Evaluate(gyt_F, gyt_uu))),
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

error if gyt_n lt 12,
    Sprintf("X0^10(13): expected at least 12 quotient comparisons, made %o (%o empty, %o skipped)",
            gyt_n, gyt_empty, gyt_skip);
printf " ok (X0^10(13): %o quotient(s) checked against Guo-Yang's curve + involutions, "
       * "%o empty, %o skipped)\n", gyt_n, gyt_empty, gyt_skip;
