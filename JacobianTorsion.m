// The zeta function of a quotient X_0^D(N)/W over a finite field, and the torsion of its Jacobian,
// from the trace formula alone -- no model of the curve is used anywhere here.
//
// Why the torsion matters: a Borcherds form with integral principal part has a divisor class of
// finite order, so when the Jacobian has no rational torsion the form is a function and Guo-Yang's
// Assumption 40 is not needed; when there is torsion, the assumption is a genuine hypothesis and the
// coordinates have to be built from e-th powers. Deciding which case one is in needed a model until
// now; it does not.

declare verbose JacobianTorsion, 2;

intrinsic ZetaNumeratorViaTrace(X::ShimuraQuot, p::RngIntElt : Check := true) -> RngUPolElt
{The numerator P(T) of the zeta function of X over F_p, from the point counts of the trace formula
over F_p up to F_(p^g): the power sums of the Frobenius eigenvalues give the first g coefficients by
Newton's identities and the functional equation gives the rest. Requires p of good reduction, that
is p not dividing D*N, and X`g assigned.

With Check, the count over F_(p^(g+1)) is computed as well and compared with what P predicts. That
comparison is the test of the result: verifying the functional equation would prove nothing, since
the upper coefficients are built from it, and a wrong genus passes that test at every prime.}
    require assigned X`g: "assign X`g";
    require IsPrime(p) and (X`D * X`N) mod p ne 0: "p must be a prime of good reduction, so p must not divide D*N";
    g := X`g;
    PT<T> := PolynomialRing(Integers());
    if g eq 0 then return PT!1; end if;
    s := [Rationals() | p^k + 1 - ComputePointsViaTrace(X, p, k) : k in [1..g]];
    // Newton: k e_k = sum_{j<k} (-1)^(k-j-1) e_j s_(k-j), with e_0 = 1
    e := [Rationals() | ];
    for k in [1..g] do
        acc := Rationals()!0;
        for j in [0..k-1] do
            ej := (j eq 0) select Rationals()!1 else e[j];
            acc +:= (-1)^(k-j-1) * ej * s[k-j];
        end for;
        Append(~e, acc/k);
    end for;
    cf := [Rationals() | 1] cat [(-1)^k * e[k] : k in [1..g]];
    for k in [g+1..2*g] do Append(~cf, p^(k-g) * cf[2*g-k+1]); end for;
    require &and[IsIntegral(c) : c in cf]:
        Sprintf("the point counts at p = %o do not come from a curve of genus %o: the zeta numerator "
                * "has non-integral coefficients %o", p, g, cf);
    P := PT![Integers() | c : c in cf];
    if Check then
        // the power sums the polynomial predicts, by Newton's identities read the other way, as far
        // as the first exponent the construction did not use
        ee := [Rationals() | (-1)^k * Coefficient(P, k) : k in [1..2*g]];
        ss := [Rationals() | ];
        for k in [1..g+1] do
            acc := (-1)^(k-1) * k * ee[k];
            for i in [1..k-1] do acc +:= (-1)^(i-1) * ee[i] * ss[k-i]; end for;
            Append(~ss, acc);
        end for;
        require ss[g+1] eq p^(g+1) + 1 - ComputePointsViaTrace(X, p, g+1):
            Sprintf("the zeta numerator at p = %o built from the counts up to F_(p^%o) does not "
                    * "predict the count over F_(p^%o): the counts do not come from a curve of "
                    * "genus %o", p, g, g+1, g);
    end if;
    return P;
end intrinsic;

intrinsic JacobianTorsionBound(X::ShimuraQuot : MaxPrime := 60, NumPrimes := 8, Check := true) -> RngIntElt, SeqEnum
{A multiple of the order of the rational torsion of the Jacobian of X, with the pairs <p, #J(F_p)>
it was computed from. For an abelian variety over Q with good reduction at an odd prime p, reduction
is injective on the rational torsion (Katz, Galois properties of torsion points on abelian varieties,
Invent. Math. 62 (1981), Appendix; Silverman, AEC VII.3.1 for elliptic curves), so the torsion order
divides every #J(F_p) = P_p(1) and hence their greatest common divisor. The counts come from the
trace formula, so no model of X is needed. A bound of 1 proves the Jacobian torsion-free, which is
what makes a Borcherds form with integral principal part a function and so removes Guo-Yang's
Assumption 40; a bound of 2 that no prime improves is evidence of a rational point of order 2, since
by Chebotarev a 2-torsion subgroup with no Galois-fixed point would give some prime an odd order.
Check is passed to ZetaNumeratorViaTrace, which then counts over one further extension per prime.}
    require assigned X`g and X`g ge 1: "assign X`g, which must be at least 1";
    orders := [];
    for p in PrimesInInterval(3, MaxPrime) do
        if (X`D * X`N) mod p eq 0 then continue; end if;
        n := Evaluate(ZetaNumeratorViaTrace(X, p : Check := Check), 1);
        Append(~orders, <p, n>);
        vprintf JacobianTorsion: "p = %o: #J(F_p) = %o\n", p, n;
        if #orders ge NumPrimes then break; end if;
    end for;
    require #orders gt 0: "no prime of good reduction below MaxPrime";
    return GCD([t[2] : t in orders]), orders;
end intrinsic;
