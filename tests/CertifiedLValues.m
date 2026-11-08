// The certified central L-values of HeegnerHeights.m against known values: L'(37a1, 1) and L(11a1, 1) to
// 35 digits (LMFDB, curves 37.a1 and 11.a2, special value of the L-function), each inside the proved error
// bound, which must be small; and the BSD leading term of 37a1 agreeing with L'(37a1, 1).

printf "Testing the certified L-values at s = 1...";
R := RealField(60);
E := EllipticCurve("37a1");
v, e := CertifiedLDerivativeAtOne(E : Precision := 50);
assert e lt 10^-24 and Abs(v - R!(30599977383405230182048368332167647445/10^38)) le e + 10^-35;
v, e := CertifiedLValueAtOne(EllipticCurve("11a1") : Precision := 50);
assert e lt 10^-24 and Abs(v - R!(25384186085591068433775892335090946105/10^38)) le e + 10^-35;
assert Abs(BSDLeadingTermRankOne(E) - R!(30599977383405230182048368332167647445/10^38)) lt 10^-20;
printf " ok (L'(37a1, 1), L(11a1, 1), BSD for 37a1)\n";
