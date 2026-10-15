// Independent prediction from Elkies (as recorded in tests/CMFieldsElkies.m): the function field of
// X_0^6(1) is Q(t)(sqrt(-t), sqrt(3(t-1))), with the star Hauptmodul t taking 0, 1, oo at the
// elliptic points of X*(6) of orders 2, 4, 6 (CM discriminants -24, -4, -3).
// Set y^2 = -t, u^2 = 3(t-1): then u^2 = -3(y^2+1); with u = i*sqrt(3)*v this is v^2 = y^2+1,
// parametrised by m via v = (m+1/m)/2, y = (m-1/m)/2, so t = -y^2 = -((m^2-1)/(2m))^2.
// The four elliptic points of X_0^6(1) have CM by Z[i] (disc -4, order 2) and Z[rho] (disc -3,
// order 3); disc -3 is t = oo and disc -4 is t = 1.
P<m> := FunctionField(Rationals());
t := -((m^2-1)/(2*m))^2;
printf "t(m) = %o, degree %o\n", t, Max(Degree(Numerator(t)), Degree(Denominator(t)));
// t = oo  <=>  m = 0 or oo   (the two order-3 points)
// t = 1   <=>  solve
rts := Roots(Numerator(t - 1));
printf "t = 1 at m = %o\n", rts;
rts0 := Roots(Numerator(t));
printf "t = 0 at m = %o  (disc -24, fixed by w_6, NOT elliptic for Gamma)\n", rts0;
// over Q(i)
K<i> := QuadraticField(-1);
PK<mm> := PolynomialRing(K);
tt := -((mm^2-1))^2 + (2*mm)^2*0;
printf "t = 1 roots over Q(i): %o\n", Roots(PK ! (-( (mm^2-1)^2 ) - (2*mm)^2));
// cross-ratio of (order-3 pair ; order-2 pair) = (0, oo ; i, -i)
CR := func< z1,z2,z3,z4 | ((z1-z3)*(z2-z4))/((z1-z4)*(z2-z3)) >;
KK<ii> := QuadraticField(-1);
printf "cross-ratio (0, oo ; i, -i) = limit = %o\n", (-ii)/(ii);
quit;
