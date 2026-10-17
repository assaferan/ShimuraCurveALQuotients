SetColumns(0);
AttachSpec("ShimuraQuotients.spec");
P<x> := PolynomialRing(Rationals());
models := AssociativeArray();
// the fitted genus-3 candidate
models[[Integers()|1,11,39,429]] := [* <3, -3*x^7-22*x^6-55*x^5-69*x^4-55*x^3-22*x^2-3*x, x^4+x^3+x^2+x+1> *];
// the fitted star (genus 1)
models[[Integers()|1,3,11,13,33,39,143,429]] := [* <1, -25447851/236421376*x^4 - 573417/476656*x^3 + 18083/7688*x^2 + 13*x + 13, P!0> *];
nchk, nfail := VerifyModelSet(models, 39, 11 : NPrimes := 6);
printf "\nVerifyModelSet(39,11 fitted): %o checks, %o failures\n", nchk, nfail;

// Jacobian of the star
Cst := HyperellipticCurve(-25447851/236421376*x^4 - 573417/476656*x^3 + 18083/7688*x^2 + 13*x + 13);
printf "star model genus = %o\n", Genus(Cst);
pts := Points(Cst : Bound := 1000);
printf "rational points found on star (bd 1000): %o\n", #pts;
if #pts gt 0 then
  E, mp := EllipticCurve(Cst, pts[1]);
  E := MinimalModel(E);
  printf "star Jacobian/elliptic model: %o\n", E;
  printf "  conductor = %o\n", Conductor(E);
  printf "  Cremona = %o\n", CremonaReference(E);
end if;
printf "DONE\n";
exit;
