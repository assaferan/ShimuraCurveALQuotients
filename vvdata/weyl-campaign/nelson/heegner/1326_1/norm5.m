SetColumns(0);
load "/private/tmp/claude-501/-Users-assaferan-Documents-GitHub-ShimuraCurveALQuotients/2ccd0b0f-dcda-423c-b902-8aad70e4879f/scratchpad/heeg/octic1326.m";
PX<X> := PolynomialRing(Rationals());
b3 := -691048320879367554399/110742857919389961800; c3 := -53690826245970963317/332228573758169885400;
ok, s := IsSquare(1 + b3*c3); r1 := (1 + s)/c3; r2 := (1 - s)/c3;
Fn := &+[Coefficient(F8, k) * (r2*X - r1)^k * (X - 1)^(8-k) : k in [0..8]];
assert &and[Coefficient(Fn, k) eq 0 : k in [1,3,5,7]];
Fn *:= LCM([Denominator(c) : c in Coefficients(Fn)]);
cont := GCD([Integers() | c : c in Coefficients(Fn)]); sq := &*[Integers() | f[1]^(f[2] div 2) : f in Factorisation(cont)]; Fn := Fn div sq^2;
G := PX![Coefficient(Fn, 2*k) : k in [0..4]];
size := func<Q | &+[#IntegerToString(Integers()!c) : c in Coefficients(Q)]>;
function scaled(G, lam)   // u -> lam u with lam a SQUARE (a Q-isomorphism x -> sqrt(lam) x), then square content
    Gs := PX![Coefficient(G, k)*lam^k : k in [0..4]];
    Gs *:= LCM([Denominator(c) : c in Coefficients(Gs)])^2;
    gg := GCD([Integers() | c : c in Coefficients(Gs)]);
    sq := &*[Integers() | f[1]^(f[2] div 2) : f in Factorisation(gg)];
    return Gs div sq^2;
end function;
ps := {2, 3, 13, 17, 67, 11287, 583019};
best := G; changed := true;
while changed do
    changed := false;
    for p in ps, e in [-8..8] do
        if e eq 0 then continue; end if;
        Gs := scaled(best, p^(2*e));
        if size(Gs) lt size(best) then best := Gs; changed := true; end if;
    end for;
end while;
printf "G (square scalings only) = %o\n", best;
F8n := PX![IsEven(k) select Coefficient(best, k div 2) else 0 : k in [0..8]];
printf "Q-isomorphic to the raw model: %o\n", IsIsomorphic(HyperellipticCurve(F8n), HyperellipticCurve(F8));
printf "discriminant: %o\n", Factorisation(Integers()!Discriminant(F8n));
G0 := PX![68, -36, -55, -36, 68];
for c in [1,-1,2,-2,3,-3,6,-6,13,-13,17,-17,26,-26,34,-34,39,-39,51,-51,78,-78,102,-102,221,-221,442,-442,663,-663,1326,-1326] do
    for lam in [1, 2, 3, 6, 13, 17, 26, 34, 39, 51, 78, 102, 221, 442, 663, 1326] do
        F0 := PX![IsEven(k) select c*Coefficient(G0, k div 2)*lam^(k div 2) else 0 : k in [0..8]];
        if IsIsomorphic(HyperellipticCurve(F0), HyperellipticCurve(F8)) then printf "true Q-form: y^2 = %o * G0(%o x^2), G0 = 68u^4 - 36u^3 - 55u^2 - 36u + 68\n", c, lam; end if;
    end for;
end for;
PrintFile("/private/tmp/claude-501/-Users-assaferan-Documents-GitHub-ShimuraCurveALQuotients/2ccd0b0f-dcda-423c-b902-8aad70e4879f/scratchpad/heeg/model1326min.m", Sprintf("P<X> := PolynomialRing(Rationals()); F := %o;", F8n) : Overwrite := true);
printf "DONE\n";
