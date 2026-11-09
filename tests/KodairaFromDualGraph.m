// The Kodaira index of a genus-one Atkin-Lehner quotient at a prime of D from the metric dual graph
// (KodairaIndexFromDualGraph), against the elliptic curves the quotients are known to be.
//
// Calibration: X_0^142(1)*, X_0^158(1)* and X_0^38(5)* are 142a1, 158a1 and 190a1 (the Jacobian is the
// curve of the unique rational newform with the right Atkin-Lehner signs, and each isogeny class has a
// single curve); their Kodaira types from Cremona's tables are I_9, I_1 / I_8, I_1 / I_11, I_1 at the two
// primes of D.
// Identification: X_0^14(17)* lies in the class 238a, whose curves 238a1 and 238a2 have types I_14, I_2
// and I_7, I_4 at 2 and 7; X_0^94(3)* lies in 282b, with 282b1 of types I_8, I_1 and 282b2 of I_4, I_2 at
// 2 and 47. The dual graph picks 238a2 and 282b1.

procedure check(D, N, W, label)
    E := EllipticCurve(label);
    for p in PrimeDivisors(D) do
        n := KodairaIndexFromDualGraph(D, N, W, p);
        ks := KodairaSymbol(E, p);
        assert ks eq KodairaSymbol("I" cat IntegerToString(n));
    end for;
end procedure;

printf "Testing the Kodaira index from the dual graph on the three determined stars...";
check(142, 1, {1,2,71,142}, "142a1");
check(158, 1, {1,2,79,158}, "158a1");
check(38, 5, {1,2,19,38,5,10,95,190}, "190a1");
printf "Done\n";
printf "Identifying X_0^14(17)* and X_0^94(3)* within their isogeny classes...";
check(14, 17, {1,2,7,14,17,34,119,238}, "238a2");
check(94, 3, {1,2,47,94,3,6,141,282}, "282b1");
// and the other curve of each class is excluded at some prime
assert exists{p : p in [2, 7] | KodairaSymbol(EllipticCurve("238a1"), p) ne KodairaSymbol("I" cat IntegerToString(KodairaIndexFromDualGraph(14, 17, {1,2,7,14,17,34,119,238}, p)))};
assert exists{p : p in [2, 47] | KodairaSymbol(EllipticCurve("282b2"), p) ne KodairaSymbol("I" cat IntegerToString(KodairaIndexFromDualGraph(94, 3, {1,2,47,94,3,6,141,282}, p)))};
printf "Done\n";
