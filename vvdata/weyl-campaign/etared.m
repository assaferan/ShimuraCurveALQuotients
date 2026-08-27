// etared: eta with the argument reduced to the fundamental domain first.
//
// WHY.  The E-table evaluates eta at g^{-1} tau0 and at (a tau0 + b)/e.  Coset
// reps are ordered by (c:d), so the first ~M of them have c = 1 with d running
// up to M-1, and there
//     Im(g^{-1} tau0) = Im(tau0) / |tau0 + d|^2  ~  1/d^2 .
// Magma's DedekindEta sums the pentagonal series in q = e(z), whose term count
// grows like 1/sqrt(Im z), so the cost grows LINEARLY IN THE COSET INDEX and
// the E-table total is QUADRATIC in M.  Measured at M = 4620, PREC 600:
//     coset     3   Im 3.8e-1   ~20 terms      ~0 ms
//     coset   500   Im 5.3e-6   ~5272 terms     3.3 ms
//     coset  1000   Im 1.3e-6   ~10562 terms   10.0 ms
//     coset  2500   Im 2.1e-7   ~26432 terms   26.7 ms
// which is exactly the observed per-coset time growth (434, 1059, 1608, 2212,
// 2822 s per 500 cosets) and hence ~64 h for the base.
//
// THE FIX.  eta is only ever needed up to the modular transformations that
// define it, so reduce first:
//     eta(w + n) = e(n/24) eta(w)          (n integral)
//     eta(-1/w)  = sqrt(-i w) eta(w)   =>  eta(w) = eta(-1/w) / sqrt(-i w)
// Applying these until |Re w| <= 1/2 and |w| >= 1 puts w in the fundamental
// domain, where Im w >= sqrt(3)/2 and the series needs ~20 terms REGARDLESS of
// the original argument.  Cost per evaluation becomes constant, so the E-table
// becomes linear in the number of cosets.
//
// It also removes the need for the huge PREC of the underflow fix: the true
// values really are as small as 1e-378, but they now arise as a product of
// MODERATE factors (roots of unity and square roots), each carrying full
// relative precision, instead of from a catastrophically cancelling series.
// Relative accuracy is what the ratios need, so PREC ~ 120 suffices again.
//
// Branch: for w in the upper half plane, -i w has POSITIVE real part, so the
// principal square root is continuous there and the choice is unambiguous.

etaRed := function(z)
    CC := Parent(z);
    ii := CC.1;
    twopii := 2*Pi(CC)*ii;
    fac := CC!1;
    w := z;
    cnt := 0;
    while true do
        n := Round(Re(w));
        if n ne 0 then
            w := w - n;
            fac *:= Exp(twopii*(CC!n)/24);     // eta(w) = e(n/24) eta(w - n)
        end if;
        if Re(w)^2 + Im(w)^2 lt 1 then
            fac /:= Sqrt(-ii*w);               // eta(w) = eta(-1/w) / sqrt(-i w)
            w := -1/w;
        else
            break;
        end if;
        cnt +:= 1;
        error if cnt gt 100000, "eta reduction did not terminate", z;
    end while;
    return fac * DedekindEta(w);
end function;
