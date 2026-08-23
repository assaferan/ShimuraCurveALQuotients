// The offline verdict for the degree-weighting experiment on X0^10(23):
// given the raw Schofer LogSums at the anchors (-40, -43) and the failing
// quadratic disc -68 (from /tmp/degw_10_23.log), test which power of 23 in the
// -68 cell makes one of the four sign-candidate minpolys admissible in the
// star field of definition K -- i.e. has a non-rational root in K.
//
// Fill in the three pairs of LogSums below from the experiment output, as
// RATIONAL-COEFFICIENT log data: each value is exp(sum_i q_i log p_i), which we
// carry as the exact rational prod p_i^{q_i} (all Schofer LogSums on this base
// are rational combinations of logs of primes with rational exponents; if an
// exponent is non-integral, carry the value as an element of a suitable
// radical extension or compare numerically at 200 digits instead).
//
//   magma -b snapcheck.m
AttachSpec("ShimuraQuotients.spec");
D := 10; N := 23;
Xstar := CreateShimuraQuot(D, N, Set(Divisors(D*N)));
Xstar`g := GenusShimuraCurveQuotient(D, N, Xstar`W); Xstar`CurveID := 0;
curves := GetQuotientsAndGenera([Xstar]);
_ := exists(star){c : c in curves | IsStarCurve(c)};

// star field of definition at d = -68
flds := FieldsOfDefinitionOfCMPointFast(star, -68 : MaxDegree := 2);
K := flds[1];
printf "star field at -68: %o\n", K;

// ---- FILL IN from /tmp/degw_10_23.log (RationalNumber of each LogSum, or the
// LogSum print itself translated to a rational) ----
// scale      = s(-43)   [anchor: disc where stilde = 0]
// scaletilde = stilde(-40) [anchor: disc where s = 0]
// ns68, nst68 = the -68 cells for s, stilde (raw, BEFORE anchor division)
scale := 23/(2^17*5);
scaletilde := 1/5;
ns68 := 23/(2^34*5^2);
nst68 := 1/(2^2*5^2);
error if scale eq 0, "fill in the four values from the experiment log first";

R<x> := PolynomialRing(Rationals());
for shift in [-2, -1, 0, 1, 2] do
    // hypothesis: the true -68 cells differ from the pipeline's by 23^(shift*mult)
    // with mult(-1) = 1, mult(-2) = 0 on this base: only the s-row shifts.
    norm_s := (ns68 * 23^shift) / scale^2;
    norm_st := nst68 / scaletilde^2;
    good := 0;
    for eps in [[1,1], [1,-1], [-1,1], [-1,-1]] do
        tr := 1 - eps[1]*norm_st + eps[2]*norm_s;
        pol := x^2 - tr*x + eps[2]*norm_s;
        rts := Roots(pol, K);
        if #rts ne 0 and not &and[rt[1] in Rationals() : rt in rts] then
            good +:= 1;
        end if;
    end for;
    printf "shift 23^%o on the s-row: %o admissible minpoly(s)\n", shift, good;
end for;
quit;
