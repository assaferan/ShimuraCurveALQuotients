SetColumns(0);
load "/private/tmp/claude-501/-Users-assaferan-Documents-GitHub-ShimuraCurveALQuotients/2ccd0b0f-dcda-423c-b902-8aad70e4879f/scratchpad/heeg/model210.m";
size := func<Q | &+[#IntegerToString(Integers()!c) : c in Coefficients(Q)]>;
function scl(F, lam)   // x -> lam x (always a Q-isomorphism), then clear denominators and square content
    G := P![Coefficient(F, k)*lam^k : k in [0..8]];
    G *:= LCM([Denominator(c) : c in Coefficients(G)])^2;
    gg := GCD([Integers() | c : c in Coefficients(G)]); sq := &*[Integers() | f[1]^(f[2] div 2) : f in Factorisation(gg)];
    return G div sq^2;
end function;
best := F; changed := true;
while changed do
    changed := false;
    for p in [2,3,5,7,13], e in [-10..10] do
        if e eq 0 then continue; end if;
        G := scl(best, p^e);
        if size(G) lt size(best) then best := G; changed := true; end if;
    end for;
end while;
// also x -> 1/x (reverse) and shifts by small rationals
best2 := best;
for sh in [a/b : a in [-6..6], b in [1,2,3,4,5,6]] do
    G := scl(Evaluate(best, X + sh), 1); if size(G) lt size(best2) then best2 := G; end if;
end for;
printf "model: y^2 = %o\n", best2;
printf "reversed: y^2 = %o\n", P!Reverse(Coefficients(best2));
printf "discriminant: %o\n", Factorisation(Integers()!Discriminant(best2));
printf "isomorphic to the raw model: %o\n", IsIsomorphic(HyperellipticCurve(best2), HyperellipticCurve(F));
