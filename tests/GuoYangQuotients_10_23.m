// tests/GuoYangQuotients_10_23.m -- oracle for the two cover keys of X_0^10(23) that
// EquationsByRebase filled, checked against Guo-Yang.
//
// ⚠ THIS IS A PARTIAL ORACLE, DELIBERATELY, and it covers exactly the keys that needed covering.
// 10_23's curve has genus 9 (degree 20) and CurveQuotient was OOM-killed on it, so the generic
// sweep in tests/GuoYangQuotientOracle.m excludes this base. Guo-Yang's w_2 is the non-diagonal
// Mobius map ((2x+1)/(x-2), -5^5 y/(x-2)^10), whose invariants are messy; w_2 and its compositions
// are NOT derived here. What IS derived are the two keys that were empty until 2026-09-09 --
// [1,5] and [1,46] -- so nothing newly produced goes unchecked.
//
// Guo-Yang publish w_2, w_5 = (-1/x, -y/x^10) and w_230 = (x,-y). The two easy ones:
//   * w_5  : y/x^5 is INVARIANT      ((-y/x^10) / (-1/x^5) = y/x^5), and t = x - 1/x is invariant.
//   * w_46 = w_5 w_230 = (-1/x, y/x^10) : now y/x^5 is ANTI-invariant, so pair it with the
//     anti-invariant x + 1/x; the product is invariant and (x+1/x)^2 = t^2 + 4.
// Both then need f/x^10 expressed in t, which works because f is ANTI-PALINDROMIC:
// a_{20-k} = (-1)^k a_k. Then f/x^10 = a_10 + sum_m a_{10+m} p_m where p_m = x^m + (-1/x)^m
// satisfies p_m = t*p_{m-1} + p_{m-2} (the roots x and -1/x have sum t and product -1).
// ⚠ Both the anti-palindromy and the resulting identity are VERIFIED below, not assumed.

gyw_Px<x> := PolynomialRing(Rationals());
gyw_f := -43*x^20 + 318*x^19 - 1071*x^18 + 3014*x^17 - 10540*x^16 + 28266*x^15 - 72217*x^14
         + 81478*x^13 - 62765*x^12 - 68732*x^11 + 18840*x^10 + 68732*x^9 - 62765*x^8
         - 81478*x^7 - 72217*x^6 - 28266*x^5 - 10540*x^4 - 3014*x^3 - 1071*x^2 - 318*x - 43;
error if Genus(HyperellipticCurve(gyw_f)) ne 9, "X0^10(23): Guo-Yang's curve should have genus 9";

gyw_K<X> := FunctionField(Rationals());
error if (-1/X^10)^2 * Evaluate(gyw_f, X) ne Evaluate(gyw_f, -1/X),
    "X0^10(23): w_5 = (-1/x, -y/x^10) does not preserve Guo-Yang's curve";
error if not &and[Coefficient(gyw_f,20-k) eq (-1)^k*Coefficient(gyw_f,k) : k in [0..20]],
    "X0^10(23): f is not anti-palindromic, so the p_m expansion below is invalid";

gyw_Pt<t> := PolynomialRing(Rationals());
gyw_p := [gyw_Pt| 2, t];
for gyw_m in [2..10] do Append(~gyw_p, t*gyw_p[gyw_m] + gyw_p[gyw_m-1]); end for;
gyw_P := Coefficient(gyw_f,10) + &+[ Coefficient(gyw_f,10+gyw_m) * gyw_p[gyw_m+1] : gyw_m in [1..10] ];
error if Evaluate(gyw_P, X - 1/X) ne Evaluate(gyw_f, X)/X^10,
    "X0^10(23): the identity f/x^10 = P(x - 1/x) FAILS -- the p_m recursion is wrong";

gyw_oracle := [*
  <[1,5],  HyperellipticCurve(gyw_P),               "(x - 1/x, y/x^5)">,
  <[1,46], HyperellipticCurve(gyw_P*(t^2+4)),       "(x - 1/x, (y/x^5)(x + 1/x))">
*];

function gyw_model_curve(e)
    if (#e ge 3) and (Type(e[3]) eq RngUPolElt) and (e[3] ne 0) then
        return HyperellipticCurve(e[2], e[3]);
    end if;
    return HyperellipticCurve(e[2]);
end function;

gyw_models := eval (Read("data/models/models_10_23.m") cat "\nreturn models;");
gyw_n := 0;
for gyw_o in gyw_oracle do
    gyw_lab, gyw_Cq, gyw_gen := Explode(gyw_o);
    gyw_ok, gyw_es := IsDefined(gyw_models, [Integers()| s : s in gyw_lab]);
    error if (not gyw_ok) or (#gyw_es eq 0),
        Sprintf("X0^10(23): key %o is missing or empty -- it was filled on 2026-09-09 and this "
                * "test exists to check it", gyw_lab);
    gyw_Cs := gyw_model_curve(gyw_es[1]);
    error if Genus(gyw_Cs) ne Genus(gyw_Cq),
        Sprintf("X0^10(23) W=%o: our genus %o vs Guo-Yang's %o -- wrong object",
                gyw_lab, Genus(gyw_Cs), Genus(gyw_Cq));
    error if not IsIsomorphic(gyw_Cs, gyw_Cq),
        Sprintf("X0^10(23) W=%o: our stored curve is NOT isomorphic to Guo-Yang's quotient by "
                * "the invariants %o", gyw_lab, gyw_gen);
    gyw_n +:= 1;
end for;

// X/w_230 is the bare x-line, hence a P^1 -- check it is genus 0 WITH a rational point.
gyw_ok230, gyw_es230 := IsDefined(gyw_models, [Integers()|1,230]);
if gyw_ok230 and #gyw_es230 gt 0 and Type(gyw_es230[1][2]) ne MonStgElt then
    gyw_C := gyw_model_curve(gyw_es230[1]);
    error if Genus(gyw_C) ne 0, "X0^10(23) W=[1,230]: expected genus 0";
    error if not HasRationalPoint(Conic(gyw_C)),
        "X0^10(23) W=[1,230]: X/w_230 is the x-line, so it must be a P^1 over Q";
    gyw_n +:= 1;
end if;

error if gyw_n lt 3,
    Sprintf("X0^10(23): expected at least 3 comparisons, made %o", gyw_n);
printf " ok (X0^10(23): %o quotient(s) checked against Guo-Yang -- the two keys filled on "
       * "2026-09-09, plus X/w_230; w_2's non-diagonal orbit is NOT covered)\n", gyw_n;
