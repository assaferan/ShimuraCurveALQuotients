// The Jacobian of the genus-1 quotient is y^2 + x*y = x^3 - 24x + 63, i.e. [1,0,0,-24,63].
// I previously typed [0,0,0,-24,63], dropping a_1, which is a different curve of conductor 276624.
SetColumns(0);
E := MinimalModel(EllipticCurve([1,0,0,-24,63]));
printf "Jacobian of X_0^39(11)/w_3 : %o, conductor %o\n", CremonaReference(E), Conductor(E);
ps := [2,3,5,7,11,13,17,19,23,29,31,37];
aE := [ TraceOfFrobenius(E,p) : p in ps ];
printf "   a_p = %o\n", aE;
// the W-fixed newform with quaternionic w_3 = +1
for orb in Newforms(CuspForms(Gamma0(429), 2)) do
    f := orb[1];
    if Level(Parent(f)) ne 429 or Degree(BaseRing(Parent(f))) ne 1 then continue; end if;
    q := qExpansion(f, 40);
    a3 := Integers()!Coefficient(q,3); a11 := Integers()!Coefficient(q,11); a13 := Integers()!Coefficient(q,13);
    w3 := a3; w11 := -a11; w39 := a3*a13; w429 := a3*a13*(-a11);
    if not (w11 eq 1 and w39 eq 1 and w429 eq 1 and w3 eq 1) then continue; end if;
    af := [ Integers() ! Coefficient(q,p) : p in ps ];
    printf "the W-fixed newform with w_3 = +1:\n   a_p = %o\n", af;
    printf "\nIDENTICAL to the model's elliptic quotient at all %o primes: %o\n", #ps, af eq aE;
end for;
print "CERT12 DONE";
