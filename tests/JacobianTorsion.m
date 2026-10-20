// The zeta numerator and the torsion bound from the trace formula (JacobianTorsion.m), against
// objects this repository does not produce.
//
// [1] OUTSIDE CHECK. Three star quotients of genus one whose Jacobian is known from the newform and
//     a single-curve isogeny class: X_0^142(1)* = 142a1, X_0^158(1)* = 158a1, X_0^38(5)* = 190a1.
//     For a genus-one curve the zeta numerator is 1 - a_p T + p T^2, so the polynomial computed from
//     the trace formula must equal the one Magma reads off the elliptic curve, coefficient by
//     coefficient. The three curves have trivial rational torsion (Cremona's tables, and the LMFDB
//     records them as 142.a1, 158.a1, 190.a1 with trivial torsion), so the bound must come out 1.
// [2] A CASE WITH TORSION. X_0^14(17)* lies in the isogeny class 238a, and both curves of that class
//     have rational torsion Z/2 (Cremona's tables). Whichever curve the star is, the bound must
//     therefore be even -- this is independent of the identification of the star with 238a2, which
//     is this repository's own result.
// [3] THE FUNCTIONAL EQUATION is checked inside the intrinsic at every prime, which is a control on
//     the point counts themselves: counts that do not come from a curve of the stated genus fail it.
// [4] REGRESSION, CODE-DERIVED. The two genus-three stars for which the bound was wanted:
//     X_0^6(379)* comes out 1 (so Assumption 40 is not needed there) and X_0^46(31)* comes out 2.

printf "Testing the zeta numerator and the torsion bound from the trace formula...";

PT<T> := PolynomialRing(Integers());

procedure check_genus_one(D, N, lab)
    W := {m : m in Divisors(D*N) | GCD(m, (D*N) div m) eq 1};
    X := CreateShimuraQuot(D, N, W);
    X`g := GenusShimuraCurveQuotient(D, N, W);
    assert X`g eq 1;
    E := EllipticCurve(lab);
    assert #TorsionSubgroup(E) eq 1;
    n := 0;
    for p in PrimesInInterval(3, 40) do
        if (D*N) mod p eq 0 then continue; end if;
        assert ZetaNumeratorViaTrace(X, p) eq PT![1, -TraceOfFrobenius(E, p), p];
        n +:= 1;
    end for;
    assert n ge 8;
    b := JacobianTorsionBound(X);
    assert b eq 1;
    printf " %o(%o good primes)", lab, n;
end procedure;

check_genus_one(142, 1, "142a1");
check_genus_one(158, 1, "158a1");
check_genus_one(38, 5, "190a1");

// [2] the class 238a has torsion Z/2 on both of its curves, so the bound must be even
W14 := {m : m in Divisors(14*17) | GCD(m, (14*17) div m) eq 1};
X14 := CreateShimuraQuot(14, 17, W14);
X14`g := GenusShimuraCurveQuotient(14, 17, W14);
assert X14`g eq 1;
assert {#TorsionSubgroup(C) : C in IsogenousCurves(EllipticCurve("238a1"))} eq {2};
b14 := JacobianTorsionBound(X14);
assert IsEven(b14);

// [4] the two genus-three stars, as regression values
for t in [<6, 379, 1>, <46, 31, 2>] do
    D := t[1]; N := t[2]; expected := t[3];
    W := {m : m in Divisors(D*N) | GCD(m, (D*N) div m) eq 1};
    X := CreateShimuraQuot(D, N, W);
    X`g := GenusShimuraCurveQuotient(D, N, W);
    assert X`g eq 3;
    b, orders := JacobianTorsionBound(X : MaxPrime := 60, NumPrimes := 6);
    assert b eq expected;
    assert #orders eq 6;
end for;

// [5] NEGATIVE CONTROL. The checks inside the intrinsic have to be able to fail. Told that the
// genus-one star X_0^142(1)* has genus two, the counts over F_p and F_(p^2) no longer come from a
// curve of that genus, and the integrality or the functional equation must refuse; a run that
// answered anyway would mean neither check constrains the counts.
Xc := CreateShimuraQuot(142, 1, {m : m in Divisors(142) | GCD(m, 142 div m) eq 1});
Xc`g := 2;                       // deliberately wrong: the star has genus one
refused := 0;
for p in [3, 5, 11, 13] do
    try
        _ := ZetaNumeratorViaTrace(Xc, p);
    catch e
        refused +:= 1;
    end try;
end for;
assert refused ge 3;

printf " done: three Jacobians identified against Cremona, one even bound, two genus-three bounds, "
    * "and a wrong genus refused at %o of 4 primes\n", refused;
