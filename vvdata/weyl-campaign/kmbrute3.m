// kappa-minus, step 3: confirm W_N(1) = (N-1)*ord_N(m) from the INDEPENDENT
// density, not from the code.  The banked closed form was read off the code
// itself (and tests/EisensteinLocalFactors.m pins it against that reading), so
// it has never been checked against anything external.
AttachSpec("ShimuraQuotients.spec");
alpha := function(p, A, m, kmax)
    prev := Rationals()!0;
    for k in [2..kmax] do
        pk := p^k; cnt := 0;
        for x1 in [0..pk-1] do for x2 in [0..pk-1] do
            num := A[1][1]*x1^2 + 2*A[1][2]*x1*x2 + A[2][2]*x2^2;
            if IsEven(num) and ((num div 2) - m) mod pk eq 0 then cnt +:= 1; end if;
        end for; end for;
        cur := Rationals()!cnt/pk;
        if k ge 3 and cur eq prev then return cur, true; end if;
        prev := cur;
    end for;
    return prev, false;
end function;
printf "scaled hyperbolic A = p*[[0,1],[1,0]] : alpha_p(m) vs (p-1)*ord_p(m)\n";
bad := 0; n := 0;
for p in [3,5] do
    A := Matrix(Integers(),2,2,[0,p,p,0]);
    kmax := p eq 3 select 6 else 4;
    for m in [1..(p eq 3 select 30 else 26)] do
        a, ok := alpha(p, A, m, kmax);
        if not ok then continue; end if;
        pred := (p-1)*Valuation(m, p);
        n +:= 1;
        if a ne pred then bad +:= 1; printf "  p=%o m=%-3o alpha=%o pred=%o MISMATCH\n", p, m, a, pred; end if;
    end for;
    printf "  p = %o : checked through m, running total %o bad of %o\n", p, bad, n;
end for;
printf "\n%o of %o agree with (p-1)*ord_p(m), from brute force alone\n", n-bad, n;
printf "DONE\n";
quit;
