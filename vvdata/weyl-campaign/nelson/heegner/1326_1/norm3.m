SetColumns(0);
load "/private/tmp/claude-501/-Users-assaferan-Documents-GitHub-ShimuraCurveALQuotients/2ccd0b0f-dcda-423c-b902-8aad70e4879f/scratchpad/heeg/octic1326.m";
PX<X> := PolynomialRing(Rationals());
b3 := -691048320879367554399/110742857919389961800; c3 := -53690826245970963317/332228573758169885400;
disc := 1 + b3*c3; ok, s := IsSquare(disc); printf "1 + bc square: %o\n", ok;
r1 := (1 + s)/c3; r2 := (1 - s)/c3;
Fn := &+[Coefficient(F8, k) * (r2*X - r1)^k * (X - 1)^(8-k) : k in [0..8]];
assert &and[Coefficient(Fn, k) eq 0 : k in [1,3,5,7]];
Fn *:= LCM([Denominator(c) : c in Coefficients(Fn)]);
G := PX![Coefficient(Fn, 2*k) : k in [0..4]];
// scaling u -> lambda u, lambda = prod p^e_p, to minimise the p-adic sizes, then remove square content
cont := GCD([Integers() | c : c in Coefficients(G)]); G := G div cont;
ps := Set(PrimeDivisors(Integers()!Coefficient(G,0))) join Set(PrimeDivisors(Integers()!Coefficient(G,4)));
lam := 1;
for p in ps do
    best := 0; bestv := 10^9;
    for e in [-12..12] do
        v := Max([Valuation(Coefficient(G, k), p) + k*e : k in [0..4] | Coefficient(G,k) ne 0]) - Min([Valuation(Coefficient(G, k), p) + k*e : k in [0..4] | Coefficient(G,k) ne 0]);
        m := Min([Valuation(Coefficient(G, k), p) + k*e : k in [0..4] | Coefficient(G,k) ne 0]);
        score := v + Abs(m);
        if score lt bestv then bestv := score; best := e; end if;
    end for;
    lam *:= p^best;
end for;
Gs := PX![Coefficient(G, k)*lam^k : k in [0..4]];
Gs *:= LCM([Denominator(c) : c in Coefficients(Gs)]);
gg := GCD([Integers() | c : c in Coefficients(Gs)]); sq := &*[Integers() | f[1]^(f[2] div 2) : f in Factorisation(gg)]; Gs := Gs div sq^2;
printf "quartic G (y^2 = G(x^2)): %o\n", Gs;
F8n := PX![IsEven(k) select Coefficient(Gs, k div 2) else 0 : k in [0..8]];
Hn := HyperellipticCurve(F8n); H := HyperellipticCurve(F8);
printf "Shioda invariants agree with the raw model: %o;  Q-isomorphic: %o\n", ShiodaInvariantsEqual(ShiodaInvariants(Hn), ShiodaInvariants(H)), IsIsomorphic(Hn, H);
printf "discriminant: %o\n", Factorisation(Integers()!Discriminant(F8n));
PrintFile("/private/tmp/claude-501/-Users-assaferan-Documents-GitHub-ShimuraCurveALQuotients/2ccd0b0f-dcda-423c-b902-8aad70e4879f/scratchpad/heeg/model1326.m", Sprintf("P<x> := PolynomialRing(Rationals()); F := %o;", F8n) : Overwrite := true);
// elliptic quotients of the normal form
A, mA := AutomorphismGroup(Hn); printf "automorphism group order %o\n", #A;
for a in A do
    if Order(a) ne 2 then continue; end if;
    phi := mA(a); eqs := DefiningPolynomials(phi); R := Parent(eqs[1]);
    if eqs[1] eq R.1 and eqs[3] eq R.3 then continue; end if;
    C, q := CurveQuotient(AutomorphismGroup(Hn, [phi]));
    if Genus(C) ne 1 then printf "quotient of genus %o\n", Genus(C); continue; end if;
    pts := Points(C : Bound := 500);
    if #pts gt 0 then E := MinimalModel(EllipticCurve(C, pts[1])); printf "quotient by x -> %o: %o\n", eqs[1]/eqs[3], CremonaReference(E);
    else J := Jacobian(GenusOneModel(C)); printf "quotient by x -> %o: Jacobian %o\n", eqs[1]/eqs[3], CremonaReference(MinimalModel(J)); end if;
end for;
printf "DONE\n";
