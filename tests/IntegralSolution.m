// tests/IntegralSolution.m
//
// The IntegralSolution parameter of BorcherdsForms is DEFAULT FALSE, so without this test no CI
// job ever executes the code path it adds.
//
// What it does: Schofer requires c_eta(m) in Z, but `Solution` returns ONE arbitrary point of
// the affine space sol + Kernel(coeffs_trunc).  Where that kernel is nonzero the same divisor
// may admit an INTEGRAL representative that Solution simply did not pick -- so a base's
// "non-integral principal parts" verdict can be an artifact of the choice rather than a property
// of the divisor.  Measured over the whole NONINTEGRAL class, 7 of 18 bases were such artifacts.
//
// X_0(33,2) is the canonical case and is cheap (~21 s per BorcherdsForms call):
//     default                     principal parts NON-integral, denominator 2
//     IntegralSolution := true    principal parts INTEGRAL,     denominator 1
// with the same number of forms either way.
//
// Both directions are asserted deliberately.  Checking only that the flag gives denominator 1
// would still pass if the flag became a no-op AND the default path silently started returning
// integral solutions -- so the default is pinned as non-integral too, which is what makes this a
// test of the FLAG rather than of the base.

printf "IntegralSolution.m: 33_2 non-integral by default, integral under the flag...";

is_D := 33; is_N := 2;
is_X := CreateShimuraQuot(is_D, is_N, Set(Divisors(is_D*is_N)));
is_X`g := GenusShimuraCurveQuotient(is_D, is_N, is_X`W);
is_X`CurveID := 0;
is_curves := GetQuotientsAndGenera([is_X]);
_ := exists(is_star){c : c in is_curves | IsStarCurve(c)};

// largest denominator appearing in any form's principal part at infinity
function is_maxden(fs)
    L := 1;
    for k in Keys(fs) do
        q := qExpansionAtoo(fs[k], 1);
        v := Valuation(q);
        if v lt 0 then
            for i in [v..-1] do
                L := LCM(L, Denominator(Rationals()!Coefficient(q, i)));
            end for;
        end if;
    end for;
    return L;
end function;

is_def := BorcherdsForms(is_star, is_curves : Prec := 100);
is_int := BorcherdsForms(is_star, is_curves : Prec := 100, IntegralSolution := true);

error if #Keys(is_def) ne #Keys(is_int),
    Sprintf("IntegralSolution: form count changed, %o -> %o", #Keys(is_def), #Keys(is_int));

is_dd := is_maxden(is_def);
is_di := is_maxden(is_int);

error if is_dd eq 1,
    Sprintf("IntegralSolution: the DEFAULT path returned integral parts at 33_2 (den %o); " *
            "this test can no longer distinguish the flag from the default", is_dd);
error if is_di ne 1,
    Sprintf("IntegralSolution := true failed to clear denominators at 33_2: %o (default was %o)",
            is_di, is_dd);

printf "Done!\n";
