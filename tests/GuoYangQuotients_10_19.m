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
  <[1,10,38,95],   HyperellipticCurve(u*gyq_G),                           "(u, xz)">,
  // ⚠ w_10 AND w_95 WERE MISSING FROM THIS LIST UNTIL 2026-09-08, AND THAT OMISSION CAUSED A
  // WRONG CONCLUSION. Both quotients are CRV pairs in Guo-Yang's coordinates -- {a^2 = uF(u),
  // c^2 = uG(u)} for w_10, {y^2 = F(u), b^2 = uG(u)} for w_95 -- so they are not plain sign
  // patterns and I skipped them. Their conic c^2 = u(5u-32) HAS the rational point (0,0), so it
  // parametrises: c = t*u gives u = 32/(5-t^2), and both collapse to hyperelliptic models.
  // Because they were absent, "matches no Guo-Yang quotient" was reported for a curve that IS
  // X/w_10 -- an incomplete comparison set treated as a complete one. See gyq_w10 below.
  <[1,95],  HyperellipticCurve(16*z^8 + 960*z^6 + 41568*z^4 - 233536*z^2 - 1520), "(u, y, xz)">,
  <[1,10],  HyperellipticCurve(-512*z^6 - 33280*z^4 - 1496576*z^2 - 9728),        "(u, xy, xz)">
*];

// the calibration, asserted rather than assumed
error if Genus(gyq_oracle[1][2]) ne 2,
    "X0^10(19): X/w_190 should have genus 2 (Guo-Yang print it); the involution composition is wrong";
error if Genus(gyq_oracle[2][2]) ne 0,
    "X0^10(19): X/w_38 should have genus 0 (Guo-Yang print it); the involution composition is wrong";

// ⚠ A STORED ENTRY MAY BE <genus, f, h>, MEANING y^2 + h*y = f -- NOT y^2 = f. Nine entries
// across seven model files carry a nonzero h. Reading only e[2] silently drops it and yields a
// DIFFERENT CURVE of the same genus, which is exactly the wrong-object mistake this repo keeps
// paying for: it cost a false "defect" report against models_87_1.m on 2026-09-08, where
// 4*f + h^2 is precisely Guo-Yang's published polynomial.
function model_curve(e)
    if (#e ge 3) and (Type(e[3]) eq RngUPolElt) and (e[3] ne 0) then
        return HyperellipticCurve(e[2], e[3]);
    end if;
    return HyperellipticCurve(e[2]);
end function;

gyq_models := eval (Read("data/models/models_10_19.m") cat "\nreturn models;");
gyq_n := 0; gyq_empty := 0;
for gyq_o in gyq_oracle do
    gyq_lab, gyq_Cq, gyq_gen := Explode(gyq_o);
    gyq_ok, gyq_es := IsDefined(gyq_models, [Integers()| t : t in gyq_lab]);
    if (not gyq_ok) or (#gyq_es eq 0) then gyq_empty +:= 1; continue; end if;
    for gyq_e in gyq_es do
        if Type(gyq_e[2]) eq MonStgElt then continue; end if;      // CRV entries: not compared here
        gyq_Cs := model_curve(gyq_e);
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
error if gyq_n lt 12,
    Sprintf("X0^10(19): expected at least 12 quotient comparisons, made %o (%o keys empty)",
            gyq_n, gyq_empty);

printf " ok (X0^10(19): %o quotient(s) checked against Guo-Yang's curve + involutions, "
       * "%o key(s) still empty)\n", gyq_n, gyq_empty;
