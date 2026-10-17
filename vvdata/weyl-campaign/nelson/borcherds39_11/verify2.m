SetColumns(0);
AttachSpec("ShimuraQuotients.spec");
P<x> := PolynomialRing(Rationals());
f := -25447851/236421376*x^4 - 573417/476656*x^3 + 18083/7688*x^2 + 13*x + 13;
Cst := HyperellipticCurve(f);
Cst := IntegralModel(Cst);
Cst := ReducedModel(Cst);
printf "star integral model: %o\n", Cst;
printf "genus %o\n", Genus(Cst);
pts := Points(Cst : Bound := 10000);
printf "#rational pts (bd 10^4) = %o : %o\n", #pts, [p : p in pts][1..Minimum(4,#pts)];
if #pts gt 0 then
  E := MinimalModel(EllipticCurve(Cst, Representative(pts)));
  printf "elliptic model: %o\n  conductor %o\n", E, Conductor(E);
  try printf "  Cremona %o\n", CremonaReference(E); catch e printf "  Cremona lookup failed\n"; end try;
  // compare with 429b1
  E429 := EllipticCurve([1,0,0,-8,-16]); // placeholder; look up below
  printf "  aInvariants %o\n", aInvariants(E);
end if;
// the genus-3 candidate: Jacobian decomposition sanity via L-poly at small primes
C3 := HyperellipticCurve(-3*x^7-22*x^6-55*x^5-69*x^4-55*x^3-22*x^2-3*x, x^4+x^3+x^2+x+1);
printf "genus-3 candidate genus = %o\n", Genus(C3);
printf "DONE\n";
exit;
