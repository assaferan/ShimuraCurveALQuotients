// Probe 3.  Can Magma's IndexFormEquation do Molina's Step 4 at the indices that actually occur?
// For each published factor of degree >= 3, rebuild the monic q(x) = a^(n-1) p(x/a) and ask
// IndexFormEquation for every generator of that index; check the true one is among them, and time it.

SetColumns(0);
P<x> := PolynomialRing(Rationals());

cases := [* <"39 a", -(7*x^4+79*x^3+311*x^2+497*x+277)>,
            <"39 b", x^4+9*x^3+29*x^2+39*x+19>,
            <"55 a", -(3*x^4-32*x^3+130*x^2-237*x+163)>,
            <"55 b", x^4-8*x^3+34*x^2-83*x+81>,
            <"35 c", 172*x^3+176*x^2+60*x+7>,
            <"51 c", 7*x^3+52*x^2+116*x+68>,
            <"57 c", x^3-19*x^2+119*x-249>,
            <"65/13", 7*x^4-3*x^3-32*x^2+25*x-5>,
            <"65/5 ", 7*x^4+81*x^3+319*x^2+508*x+268>,
            <"69   ", 4*x^4-16*x^3+11*x^2+10*x+3>,
            <"85/5 ", x^4-23*x^3+183*x^2-556*x+412>,
            <"85/85", x^4+x^3-15*x^2+20*x-8> *];

for c in cases do
    nm := c[1]; p := c[2]; n := Degree(p); a := LeadingCoefficient(p);
    q := a^(n-1) * Evaluate(p, x/a);
    assert IsMonic(q) and &and[IsIntegral(co) : co in Coefficients(q)];
    M := NumberField(q);
    O := MaximalOrder(M);
    I := Integers()!Sqrt(Discriminant(q)/Discriminant(O));
    t := Cputime();
    sols := IndexFormEquation(O, I);
    el := Cputime(t);
    // is the true generator among them, up to sign and integer translation?
    gamma := M.1;
    hit := exists{s : s in sols | IsCoercible(Integers(), M!s - gamma) or IsCoercible(Integers(), M!s + gamma)};
    printf "%o: deg %o, a = %o, index = %o -> %o solutions in %.2o s, true one found: %o\n",
        nm, n, a, I, #sols, el, hit;
end for;

printf "PROBE3 DONE\n";
quit;
