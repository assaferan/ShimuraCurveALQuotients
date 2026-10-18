// Compare the two routes to X_0^210(13)/W and normalise the model by the Klein-four action (Q-isomorphisms only).
AttachSpec("ShimuraQuotients.spec"); SetColumns(0);
load "/private/tmp/claude-501/-Users-assaferan-Documents-GitHub-ShimuraCurveALQuotients/2ccd0b0f-dcda-423c-b902-8aad70e4879f/scratchpad/heeg/octic210_w15.m"; F15 := F8;
load "/private/tmp/claude-501/-Users-assaferan-Documents-GitHub-ShimuraCurveALQuotients/2ccd0b0f-dcda-423c-b902-8aad70e4879f/scratchpad/heeg/octic210_w5.m"; F5 := F8;
H15 := HyperellipticCurve(F15); H5 := HyperellipticCurve(F5);
printf "two routes Q-isomorphic: %o   (Qbar: %o)\n", IsIsomorphic(H15, H5), ShiodaInvariantsEqual(ShiodaInvariants(H15), ShiodaInvariants(H5));
H := H5; F := F5;
PX<X> := PolynomialRing(Rationals());
A, mA := AutomorphismGroup(H);
size := func<Q | &+[#IntegerToString(Integers()!c) : c in Coefficients(Q)]>;
function sqscale(G, ps)   // u -> p^(2e) u and square content: Q-isomorphisms of y^2 = G(x^2)
    best := G; changed := true;
    while changed do
        changed := false;
        for p in ps, e in [-8..8] do
            if e eq 0 then continue; end if;
            Gs := PX![Coefficient(best, k)*p^(2*e*k) : k in [0..4]];
            Gs *:= LCM([Denominator(c) : c in Coefficients(Gs)])^2;
            gg := GCD([Integers() | c : c in Coefficients(Gs)]); sq := &*[Integers() | f[1]^(f[2] div 2) : f in Factorisation(gg)]; Gs := Gs div sq^2;
            if size(Gs) lt size(best) then best := Gs; changed := true; end if;
        end for;
    end while;
    return best;
end function;
bestF := F;
for a in A do
    if Order(a) ne 2 then continue; end if;
    eqs := DefiningPolynomials(mA(a)); R := Parent(eqs[1]);
    a1 := MonomialCoefficient(eqs[1], R.1); b1 := MonomialCoefficient(eqs[1], R.3); c1 := MonomialCoefficient(eqs[3], R.1); d1 := MonomialCoefficient(eqs[3], R.3);
    if [a1,b1,c1,d1] eq [1,0,0,1] then continue; end if;
    printf "involution matrix (%o,%o,%o,%o)\n", a1, b1, c1, d1;
    if c1 eq 0 then
        // affine involution x -> -x + beta: centre at beta/2
        if a1/d1 ne -1 then continue; end if;
        beta := b1/d1; Fn := Evaluate(F, X + beta/2);
    else
        disc := (d1 - a1)^2 + 4*c1*b1; ok, s := IsSquare(disc);
        if not ok then printf "   fixed points not rational\n"; continue; end if;
        r1 := (a1 - d1 + s)/(2*c1); r2 := (a1 - d1 - s)/(2*c1);
        Fn := &+[Coefficient(F, k) * (r2*X - r1)^k * (X - 1)^(8-k) : k in [0..8]];
    end if;
    if not &and[Coefficient(Fn, k) eq 0 : k in [1,3,5,7]] then continue; end if;
    Fn *:= LCM([Denominator(c) : c in Coefficients(Fn)])^2;
    cont := GCD([Integers() | c : c in Coefficients(Fn)]); sq := &*[Integers() | f[1]^(f[2] div 2) : f in Factorisation(cont)]; Fn := Fn div sq^2;
    G := PX![Coefficient(Fn, 2*k) : k in [0..4]];
    ps := Set(PrimeDivisors(Integers()!Discriminant(Fn))) join {2,3,5,7,13};
    G := sqscale(G, ps);
    Fe := PX![IsEven(k) select Coefficient(G, k div 2) else 0 : k in [0..8]];
    assert IsIsomorphic(HyperellipticCurve(Fe), H);
    printf "even form from involution (%o,%o,%o,%o): y^2 = %o\n", a1, b1, c1, d1, Fe;
    if size(Fe) lt size(bestF) then bestF := Fe; end if;
end for;
printf "SMALLEST: y^2 = %o\n", bestF;
printf "discriminant: %o\n", Factorisation(Integers()!Discriminant(bestF));
X := CreateShimuraQuot(210, 13, {Integers() | 1,2,21,42,65,130,1365,2730}); X`g := 3;
Hb := HyperellipticCurve(bestF); bad := 0; n := 0;
for p in PrimesInInterval(11, 600) do if 2730 mod p eq 0 then continue; end if; n +:= 1; if #Points(ChangeRing(Hb, GF(p))) ne ComputePointsViaTrace(X, p, 1) then bad +:= 1; end if; end for;
printf "counts vs trace formula: %o primes, %o mismatches\n", n, bad;
for a in A do
    if Order(a) ne 2 then continue; end if;
    eqs := DefiningPolynomials(mA(a)); R := Parent(eqs[1]);
    if eqs[1] eq R.1 and eqs[3] eq R.3 then continue; end if;
    C := CurveQuotient(AutomorphismGroup(H, [mA(a)]));
    if Genus(C) eq 1 then printf "elliptic quotient: %o\n", CremonaReference(MinimalModel(Jacobian(GenusOneModel(C)))); end if;
end for;
PrintFile("/private/tmp/claude-501/-Users-assaferan-Documents-GitHub-ShimuraCurveALQuotients/2ccd0b0f-dcda-423c-b902-8aad70e4879f/scratchpad/heeg/model210.m", Sprintf("P<x> := PolynomialRing(Rationals()); F := %o;", bestF) : Overwrite := true);
printf "DONE\n";
