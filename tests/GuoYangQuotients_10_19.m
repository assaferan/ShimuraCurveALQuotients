// tests/GuoYangQuotients_10_19.m -- a COMPLETE oracle for X_0^10(19).
//
// ⚠ THE IDEA, AND IT GENERALISES TO ANY BASE WHERE GUO-YANG PRINT THE TOP CURVE. They publish only
// three quotients for 10_19 (X/<w2,w95>, X/w_190, X/w_38), so only three of our fifteen cover keys
// looked checkable. But they also publish the full curve AND the involutions -- and every quotient
// is computable from those. That turns a partial oracle into a complete one: 11 of our entries are
// checked here, against a source that is entirely external to the pipeline.
//
// Guo-Yang, Example 37:  y^2 = -8x^6+57x^4-40x^2+16,  z^2 = 5x^2-32, with
//     w_2 : (x,y,z) -> (-x, y, z),  w_5 : (x,y,z) -> (x,-y,-z),  w_19 : (x,y,z) -> (-x,-y, z).
// The action is DIAGONAL, so each quotient is the field of invariant monomials, computed below.
// Both defining polynomials are EVEN in x, so with u = x^2 they are F(u) and G(u).
//
// ⚠ CALIBRATED BEFORE BEING TRUSTED. The group law was checked against the two quotients Guo-Yang
// state in words: w_190 negates only z, leaving y^2 = f(x) of genus 2 -- which is exactly the
// equation they print for X/w_190; and w_38 negates only y, leaving z^2 = g(x) of genus 0, again
// what they print for X/w_38. A sign error in the composition would break both.

gyq_P<u> := PolynomialRing(Rationals());
gyq_F := -8*u^3 + 57*u^2 - 40*u + 16;
gyq_G := 5*u - 32;
gyq_Pz<z> := PolynomialRing(Rationals());
gyq_f := Evaluate(gyq_F, z^2);
gyq_g := Evaluate(gyq_G, z^2);
gyq_sub := (z^2 + 32)/5;                       // u in terms of z, from z^2 = G(u)

// <cover key, quotient curve, the invariant generators it comes from>
gyq_oracle := [*
  <[1,190],        HyperellipticCurve(gyq_f),                             "(x, y)">,
  <[1,38],         HyperellipticCurve(gyq_g),                             "(x, z)">,
  <[1,5],          HyperellipticCurve(gyq_f*gyq_g),                       "(x, yz)">,
  <[1,2],          HyperellipticCurve(Evaluate(gyq_F, gyq_sub)),          "(z, y)">,
  <[1,19],         HyperellipticCurve(gyq_sub*Evaluate(gyq_F, gyq_sub)),  "(z, xy)">,
  <[1,2,5,10],     HyperellipticCurve(gyq_F*gyq_G),                       "(u, yz)">,
  <[1,2,19,38],    HyperellipticCurve(gyq_G),                             "(u, z)">,
  <[1,2,95,190],   HyperellipticCurve(gyq_F),                             "(u, y)">,
  <[1,5,19,95],    HyperellipticCurve(u*gyq_F*gyq_G),                     "(u, xyz)">,
  <[1,10,19,190],  HyperellipticCurve(u*gyq_F),                           "(u, xy)">,
  <[1,10,38,95],   HyperellipticCurve(u*gyq_G),                           "(u, xz)">
*];

// the calibration, asserted rather than assumed
error if Genus(gyq_oracle[1][2]) ne 2,
    "X0^10(19): X/w_190 should have genus 2 (Guo-Yang print it); the involution composition is wrong";
error if Genus(gyq_oracle[2][2]) ne 0,
    "X0^10(19): X/w_38 should have genus 0 (Guo-Yang print it); the involution composition is wrong";

gyq_models := eval (Read("data/models/models_10_19.m") cat "\nreturn models;");
gyq_n := 0; gyq_empty := 0;
for gyq_o in gyq_oracle do
    gyq_lab, gyq_Cq, gyq_gen := Explode(gyq_o);
    gyq_ok, gyq_es := IsDefined(gyq_models, [Integers()| t : t in gyq_lab]);
    if (not gyq_ok) or (#gyq_es eq 0) then gyq_empty +:= 1; continue; end if;
    for gyq_e in gyq_es do
        if Type(gyq_e[2]) eq MonStgElt then continue; end if;      // CRV entries: not compared here
        gyq_Cs := HyperellipticCurve(gyq_e[2]);
        error if Genus(gyq_Cs) ne Genus(gyq_Cq),
            Sprintf("X0^10(19) W=%o: our genus %o vs Guo-Yang's %o -- wrong object",
                    gyq_lab, Genus(gyq_Cs), Genus(gyq_Cq));
        gyq_isok := false;
        if Genus(gyq_Cq) eq 0 then
            gyq_isok := IsIsomorphic(Conic(gyq_Cs), Conic(gyq_Cq));
        else
            gyq_isok := IsIsomorphic(gyq_Cs, gyq_Cq);
        end if;
        error if not gyq_isok,
            Sprintf("X0^10(19) W=%o: our stored curve is NOT isomorphic to Guo-Yang's quotient "
                    * "by the invariants %o", gyq_lab, gyq_gen);
        gyq_n +:= 1;
    end for;
end for;

// ⚠ COUNT THE COMPARISONS. A key that stops being produced must not turn this green silently.
error if gyq_n lt 11,
    Sprintf("X0^10(19): expected at least 11 quotient comparisons, made %o (%o keys empty)",
            gyq_n, gyq_empty);

printf " ok (X0^10(19): %o quotient(s) checked against Guo-Yang's curve + involutions, "
       * "%o key(s) still empty)\n", gyq_n, gyq_empty;
